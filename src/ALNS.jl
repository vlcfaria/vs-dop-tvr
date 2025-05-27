module ALNS

import IterTools as itr
using Random

include("AcceleratedDubins.jl")
include("Helper.jl")

struct ALNSSolution
    seq::Vector{Tuple{Int64,Int64,Int64}}
    nodes_out::Set{Int64}

    #Random solution
    function ALNSSolution(op)
        speeds, headings = op.graph.num_speeds, op.graph.num_headings
        starting_depot = rand(op.depots)

        #Get nodes and remove chosen depot
        nodes = collect(1:length(op.coordinates))
        deleteat!(nodes, starting_depot) #Remove depot
        shuffle!(nodes)

        seq = [(starting_depot, rand(1:speeds), rand(1:headings))] #With starting depot
        nodes_out = Set()
        time = 0
        while length(nodes) > 0 #Keep adding till we cant
            nxt = pop!(nodes)
            tup = (nxt, rand(1:speeds), rand(1:headings))
            new_time = time + Helper.get_dist(op.graph.graph, seq[end], tup) #Go to tup
            if new_time + Helper.get_dist_to_depot(op, tup, starting_depot) <= op.tmax #Need to go back
                push!(seq, tup) #All ok
                time = new_time
            else
                push!(nodes_out, nxt) #Take out of solution
            end 
        end
        
        return new(seq, nodes_out)
    end
end

#Returns (legal,score)
function check_legal_solution(op, sol, ini_time = 0, ini_score = 0, idx = 2)
    score, time = ini_score, ini_time
    graph = op.graph.graph

    #Check if solution doesnt start invalid
    if idx > 1 && time + Helper.get_dist_to_depot(op, sol.seq[idx-1], sol.seq[1][1]) > op.tmax
        return false, -Inf
    end

    for i in idx:length(sol.seq)
        #Resolve incoming edge to idx
        time += Helper.get_dist(graph, sol.seq[i-1], sol.seq[i])

        #Check if visiting idx is legal
        if time + Helper.get_dist_to_depot(op, sol.seq[i], sol.seq[1][1]) > op.tmax
            return false, -Inf
        end

        #Add to score
        score += op.functions[sol.seq[i][1]](time)
    end

    return true, score
end

# Destroy operators ---

#Iteratively removes a randomly selected node from a solution
function random_removal(op, sol, num_removals::Int64,)
    for _ in 1:num_removals
        #Pick a random index of the solution
        i = rand(2:length(sol.seq))
        push!(sol.nodes_out, sol.seq[i][1])
        deleteat!(sol.seq, i)
    end
end

#Iteratively removes a node, sampling by the ones that contribute less to the solution
#Higher values of p = more deterministic

#TODO make so this doesnt actually take into account only the score of the individual node, but the total node score!
#since scores are time sensitive, we want to actually figure out the node that removing the node wields the highest score
function worst_removal(op, sol, num_removals::Int64, p::Int64)
    graph = op.graph.graph

    for _ in 1:min(num_removals, length(sol.seq)-1)
        len = length(sol.seq)
        scores = Array{Tuple{Float64,Int64}}(undef, len-1) #TOTAL score of the solution when removing node i+1
        
        running_time = 0
        running_score = 0
        for i in 2:(len-1) #Will try removing i
            #Skip i: (i-1) -> (i+1)
            time = running_time + Helper.get_dist(graph, sol.seq[i-1], sol.seq[i+1])
            
            #Score seq[i+1] and continue
            score = op.functions[sol.seq[i+1][1]](time)
            _, score = check_legal_solution(op, sol, time, score, i+2)

            scores[i-1] = (score,i)
            
            #Continue by not skipping i
            running_time += Helper.get_dist(graph, sol.seq[i-1], sol.seq[i])
            running_score += op.functions[sol.seq[i][1]](running_time)
        end

        #Running time & score contain final solution skipping the last node
        scores[len-1] = (running_score, len)

        sort!(scores, rev=true) #TODO this can be optimized into a linear selection algorithm
        y = rand()
        #Pick the index, bias towards points that still wield high scores when ignored
        selected = floor(Int64, y^p * (len-1)) + 1

        _, selected_pos = scores[selected]
        selected_node = sol.seq[selected_pos][1]

        println((selected_node, scores, sol.seq))

        #Delete from the solution
        push!(sol.nodes_out, selected_node)
        deleteat!(sol.seq, selected_pos)
    end

end

#Picks n waypoints and shuffles them. Randomly remove nodes if new solution is not valid
function random_waypoint_shuffle(op, sol, num_shuffles::Int64)
    to_shuffle = rand(1:length(sol.seq), num_shuffles) #TODO potential repeats here? Maybe try without replacements

    #Shuffle!
    for i in to_shuffle
        sol.seq[i] = (sol.seq[i][1], rand(1:op.graph.num_speeds), rand(1:op.graph.num_headings))
    end

    #While violates tmax, remove randomly
    valid, _ = check_legal_solution(op,sol)
    while !valid
        victim = rand(1:length(sol.seq))

        push!(sol.nodes_out, sol.seq[victim][1])
        deleteat!(sol.seq, victim)
        valid, _ = check_legal_solution(op,sol)
    end
end

#Repair operators ---

#Iteratively inserts nodes which maximizes the total cost, stops when no legal move remaining
function greedy_insertion(op, sol)
    graph = op.graph.graph
    while length(sol.nodes_out) > 0 #Still nodes to add
        best_insertion, best_point, best_score = (-1,-1,-1), -1, -Inf

        #Iterate through solution
        running_score = 0
        running_time = 0
        for i_point in 2:(length(sol.seq)+1) #Also try including on last position
            #Try inserting all points, with all v/h
            for (node, v, h) in itr.product(sol.nodes_out,1:op.graph.num_speeds,1:op.graph.num_headings)
                #Resolve incoming edge to node
                time = running_time + Helper.get_dist(graph, sol.seq[i_point-1], (node,v,h))

                #Check legality of inserting
                if time + Helper.get_dist_to_depot(op, (node,v,h), sol.seq[1][1]) > op.tmax
                    continue
                end
                score = running_score + op.functions[node](time)
                
                #Resolve outcoming edge and test legality of our solution (the FULL sequence must abide to TMAX!)
                if i_point <= length(sol.seq)
                    time += Helper.get_dist(graph, (node,v,h), sol.seq[i_point])
                    score += op.functions[sol.seq[i_point][1]](time)
                    valid, score = check_legal_solution(op, sol, time, score, i_point+1)
                else #We stop evaluating here (no outcoming edge), so this is valid due to earlier check
                    valid = true
                end

                if valid && score > best_score
                    best_score = score
                    best_point = i_point
                    best_insertion = (node,v,h)
                end
            
                #Increase running score/time, if not the last
                if i_point <= length(sol.seq)
                    running_time += Helper.get_dist(graph, sol.seq[i_point-1], sol.seq[i_point])
                    running_score += op.functions[sol.seq[i_point][1]](running_time)
                end
            end
        end
        
        #Cant extend/insert
        if best_score == -Inf
            break
        end

        #Insert
        insert!(sol.seq, best_point, best_insertion)
        delete!(sol.nodes_out, best_insertion[1])
    end
end

#Focuses on elements for which the choice of insertion position is critical. 
#High k-regret = if best is not chosen, others will be worse
function k_regret_insertion(op,sol, k=2)
    #Calculate solution with biggest regret position, or the one with lowest possible insertions
    #One insertion is one position to insert, so ignore (speed,heading)
end

#Randomly inserts nodes. Randomization must be split -> determine which node/position will be inserted, then which heading/speed
function random_insertion(op,sol)
end

#Greedily inserts nodes that maximize the total efficiency (CLOSE TO SHAW REMOVAL?)
function greedy_efficient_insertion(op,sol)
end

end