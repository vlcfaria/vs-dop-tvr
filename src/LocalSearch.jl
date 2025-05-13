module LocalSearch #Deterministic local search procedures

import IterTools as itr

include("AcceleratedDubins.jl")
include("Helper.jl")

function get_dist(graph, p1, p2)
    return graph[p1[1], p2[1], p1[2], p2[2], p1[3], p2[3]]
end

function get_dist_to_depot(op, p1, depot)
    return op.dists_to_depot[p1[1], p1[2], p1[3], depot]
end

function score_from_running(seq::Vector{Tuple{Int64,Int64,Int64}}, op, r_score, r_time, start_idx)
    #Calculates the score, starting from the start index, and considering the r_score and r_time as the inital score and time
    elapsed_time = r_time
    score = r_score

    for i in start_idx:length(seq)
        prev = i-1

        #Cost from going from prev to i
        new_cost = op.graph.graph[seq[prev][1], seq[i][1], seq[prev][2], seq[i][2], seq[prev][3], seq[i][3]]

        #If cannot go back to depot in time, next target cannot be added, sequence is over
        if op.dists_to_depot[seq[i][1], seq[i][2], seq[i][3], 1] + new_cost + elapsed_time > op.tmax
            # PREV goes back to depot
            elapsed_time += op.dists_to_depot[seq[prev][1], seq[prev][2], seq[prev][3], 1]
            return score #i is the first index not in solution
        end

        elapsed_time += new_cost
        score += op.functions[seq[i][1]](elapsed_time)
    end

    #Got to the end of vector, go back to depot anyway
    elapsed_time += op.dists_to_depot[seq[end][1], seq[end][2], seq[end][3], 1]

    return score
end

function waypoint_change(seq::Vector{Tuple{Int64,Int64,Int64}}, op, limit_idx::Int64)
    best_change = (-1,-1,-1)
    best_pos = -1
    best_score = -Inf
    graph = op.graph

    #Special case -> Analyze the first index of the sequence
    og = seq[1]
    for (v,h) in itr.product(1:graph.num_speeds, 1:graph.num_headings)
        seq[1] = (og[1], v, h)
        score = score_from_running(seq, op, 0, 0, 2)

        if score > best_score
            best_change = seq[1]
            best_score = score
            best_pos = 1
        end
    end
    seq[1] = og

    running_score = 0 #Score obtained by until previously analyzed
    running_time = 0 #Time to reach previously analyzed

    for n in 2:limit_idx #Changing more than limit_idx wont affect
        og = seq[n]
        #Test all possibilites + calculating score
        for (v,h) in itr.product(1:graph.num_speeds, 1:graph.num_headings)
            seq[n] = (og[1], v, h)
            #Start by evaluating this
            score = score_from_running(seq, op, running_score, running_time, n)

            if score > best_score
                best_change = seq[n]
                best_score = score
                best_pos = n
            end
        end

        #Add to running score, revert changes
        seq[n] = og
        running_time += get_dist(graph.graph, seq[n-1], seq[n])
        running_score += op.functions[seq[n][1]](running_time)
    end

    #Apply best move
    seq[best_pos] = best_change

    return best_score
end

function one_point_move(seq::Vector{Tuple{Int64,Int64,Int64}}, op, limit_idx::Int64)
    best_config = (-1,-1,-1)
    best_og_pos = -1
    best_dest_pos = -1
    best_score = -Inf

    graph = op.graph
    len = length(seq)

    for og_pos in 2:len #depot cannot be changed
        og = seq[og_pos]
        deleteat!(seq, og_pos)

        #Need to get new limit idx, if the one we extracted was up to limit idx
        if og_pos <= limit_idx
            _, _, limit = Helper.calculate_seq_results(op, seq)
        else #limit is unchanged, since "real" sequence wasnt changed
            limit = limit_idx
        end

        #Running is unique to this
        running_score = 0
        running_time = 0
        for i in 2:limit #Try inserting into every position
            #Try every combination
            for (v,h) in itr.product(1:graph.num_speeds, 1:graph.num_headings)
                #Incoming edge to this
                time = running_time + get_dist(graph.graph, seq[i-1], (og[1], v, h))
                if time + get_dist_to_depot(op, (og[1], v, h), seq[1][1]) > op.tmax
                    continue #Illegal
                end
                score = running_score + op.functions[og[1]](time)

                #Outcoming edge
                time += get_dist(graph.graph, (og[1], v, h), seq[i])
                if time + get_dist_to_depot(op, seq[i], seq[1][1]) <= op.tmax #Can be extended
                    score += op.functions[seq[i][1]](time)
                    score = score_from_running(seq, op, score, time, i+1)
                end

                if score > best_score
                    best_config = (og[1], v, h)
                    best_og_pos = og_pos
                    best_dest_pos = i
                    best_score = score
                end
            end
            
            running_time += get_dist(graph.graph, seq[i-1], seq[i])
            running_score += op.functions[seq[i][1]](running_time)
        end

        #Put back
        insert!(seq, og_pos, og)
    end

    #Apply best move
    deleteat!(seq, best_og_pos)
    insert!(seq, best_dest_pos, best_config)

    return best_score
end

function one_point_exchange(seq::Vector{Tuple{Int64,Int64,Int64}}, op, limit_idx::Int64)
    best_p1_config = (-1,-1,-1)
    best_p2_config = (-1,-1,-1)
    best_p1 = -1
    best_p2 = -1
    best_score = -Inf
    graph = op.graph
    len = length(seq)

    #Find best move
    running_score = 0
    running_time = 0

    for p1 in 2:limit_idx #Changing more than limit has no effect
        og_p1 = seq[p1]
        for p2 in p1:len #Also try changing with itself
            og_p2 = seq[p2]
            #Try all pairs
            #TODO possible optimization -> if new p1 is over the new limit, v1 and h1 and irrelevant (?)
            for (v1, v2, h1, h2) in itr.product(graph.num_speeds, graph.num_speeds, 
                                                graph.num_headings, graph.num_headings)
                seq[p1] = (og_p2[1], v2, h2)
                seq[p2] = (og_p1[1], v1, h1)
                score = score_from_running(seq, op, running_score, running_time, p1)

                if score > best_score
                    best_p1 = p1
                    best_p2 = p2
                    best_p1_config = (og_p1[1], v1, h1)
                    best_p2_config = (og_p2[1], v2, h2)
                    best_score = score
                end
            end

            #Return p2 to original
            seq[p2] = og_p2
        end

        #Return p1
        seq[p1] = og_p1
        running_time += get_dist(graph.graph, seq[p1-1], seq[p1])
        running_score += op.functions[seq[p1][1]](running_time)
    end

    #Apply best move
    seq[best_p1] = best_p2_config
    seq[best_p2] = best_p1_config

    return best_score
end

#limit_idx is the first index of the solution that is not visited, or the last index, in case all points is visited
function search(sequence::Vector{Tuple{Int64,Int64,Int64}}, op, l::Int64)
    if l == 1
        method = waypoint_change
    elseif l == 2
        method = one_point_move
    else
        method = one_point_exchange
    end

    #Call method repeatedly until no better solution is found
    score, time, limit = Helper.calculate_seq_results(op, sequence)
    while true
        if score == method(sequence, op, limit)
            break
        end

        #Get new score, time and limit
        score, time, limit = Helper.calculate_seq_results(op, sequence)
    end

    return sequence, score, time
end

end