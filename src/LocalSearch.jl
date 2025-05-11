module LocalSearch #Deterministic local search procedures

import IterTools as itr

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
    len = length(seq)

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

    return seq
end

function one_point_move(sequence::Vector{Tuple{Int64,Int64,Int64}}, op, limit_idx::Int64)
    graph = op.graph
    len = length(sequence)
    if len <= 2 #Cant apply operator since first index cannot be moved
        return sequence, -1
    end

    idx = rand(1:len)
    val = sequence[idx]
    deleteat!(sequence, idx)

    if idx == 1 #If changing first point (depot), keep it there
        new_idx = 1
    elseif idx <= limit_idx #Any chosen index will modify the sequence
        new_idx = rand(2:len)
        #Assure index is different
        while new_idx == idx
            new_idx = rand(2:len)
        end
    else #Chosen index is out of the limit index, new index must be up to limit
        new_idx = rand(2:limit_idx)
    end

    #Randomly change waypoint too
    val = (val[1], rand(1:graph.num_speeds), rand(1:graph.num_headings))
    insert!(sequence, new_idx, val)

    return sequence
end

function one_point_exchange(sequence::Vector{Tuple{Int64,Int64,Int64}}, op, limit_idx::Int64)
    graph = op.graph
    len = length(sequence)
    if len <= 2 #Cant apply operator since first index cannot be moved
        return sequence, -1
    end

    idx1 = rand(2:len)

    if idx1 <= limit_idx
        idx2 = rand(2:len) #Choosing anything will modify the sequence
        while idx1 == idx2 #While also assuring it is different
            idx2 = rand(2:len)
        end
    else
        #Need to choose a value within limit index
        idx2 = rand(2:limit_idx)
    end

    #Swap and assign random waypoints
    sequence[idx1], sequence[idx2] = sequence[idx2], sequence[idx1]

    sequence[idx1] = (sequence[idx1][1], rand(1:graph.num_speeds), rand(1:graph.num_headings))
    sequence[idx2] = (sequence[idx2][1], rand(1:graph.num_speeds), rand(1:graph.num_headings))

    return sequence
end

#limit_idx is the first index of the solution that is not visited, or the last index, in case all points is visited
function search(sequence::Vector{Tuple{Int64,Int64,Int64}}, op, l::Int64, limit_idx::Int64)
    if l == 1
        return waypoint_change(sequence, op, limit_idx)
    elseif l == 2
        return one_point_move(sequence, op, limit_idx)
    else
        return one_point_exchange(sequence, op, limit_idx)
    end
end

end