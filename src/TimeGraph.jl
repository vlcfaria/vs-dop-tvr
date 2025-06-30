module TimeGraph

import IterTools as itr
using DataStructures

include("AcceleratedDubins.jl")
include("Helper.jl")
include("VsDopTvr.jl")

function build_graph_rounded(instance_path::String, time_step = 1., vehicle_params = nothing, graph_params = nothing)
    points, scores, depots, tmax = Helper.read_optvr_file(instance_path)
    
    if vehicle_params === nothing
        vehicle_params = VsDopTvr.get_cessna172_params()
    end

    if graph_params === nothing
        graph_params = VsDopTvr.DOPGraph(length(points), vehicle_params)
    end

    #Build graph
    graph_params = Helper.compute_trajectories(points, graph_params)

    #Manually traverse graph rounding values up to the ceil timestep
    for idx in eachindex(graph_params.graph)
        graph_params.graph[idx] = time_step * ceil(graph_params.graph[idx] / time_step)
    end

    return VsDopTvr.OpParameters(graph_params, points, scores, depots, tmax, Helper.compute_fastests_paths_to_depot(graph_params, depots))
end

#Topologically sorts the vertexes
function topological_sort(op)
    sorted = Vector{Tuple{Float64, Tuple{Int64,Int64,Int64}}}()
    heap = BinaryMinHeap{Tuple{Float64, Tuple{Int64, Int64, Int64}}}()
    in_heap = Set{Tuple{Float64, Tuple{Int64, Int64, Int64}}}()

    #TODO ALLOW MULTIPLE DEPOTS
    starting_depot = op.depots[1]

    #Only the depots are at t=0
    for d in op.depots
        for (v,h) in itr.product(1:op.graph.num_speeds, 1:op.graph.num_headings)
            push!(heap, (0., (d,v,h)))
            push!(in_heap, (0., (d,v,h)))
        end
    end

    #Process time wise -> This is the same as topological sort (?) Since travel time > 0
    while !isempty(heap)
        time, node = pop!(heap)
        push!(sorted, (time, node))

        #Check all neighbors
        for dest_node in itr.product(1:length(op.coordinates), 1:op.graph.num_speeds, 1:op.graph.num_headings)
            if dest_node[1] == node[1] || dest_node[1] == starting_depot
                continue
            end

            cost = Helper.get_dist(op.graph.graph, node, dest_node)

            if ((time+cost, dest_node) in in_heap) #Already in heap
                continue
            end

            if (time + cost + Helper.get_dist_to_depot(op, dest_node, starting_depot) <= op.tmax)

                push!(in_heap, (time+cost, dest_node))
                push!(heap, (time+cost, dest_node))
            end
        end
    end

    return sorted
end

function maximal_profit_path(op, sorted)
    node_sum = Dict([(val, 0.) for val in sorted]) #Track sums
    starting_depot = op.depots[1]

    #Track paths
    node_parent = Dict([(val, (-1., (-1,-1,-1))) for val in sorted])

    for (time, node) in sorted
        #Check out sucessors
        for target in 1:length(op.coordinates)
            #Check if we can go from node -> target, targetwise
            if target == starting_depot || target == node[1]
                continue
            end

            #Now, check if this node is not already visited
            parent = node_parent[(time,node)]
            while parent[1] != -1 && parent[2][1] != target
                parent = node_parent[parent]
            end

            if parent[1] != -1 
                continue 
            end

            for (v,h) in itr.product(1:op.graph.num_speeds, 1:op.graph.num_headings)
                cost = Helper.get_dist(op.graph.graph, node, (target,v,h))
                #Valid neighbor check
                if (time + cost + Helper.get_dist_to_depot(op, (target,v,h), starting_depot) > op.tmax)
                    continue
                end
                
                dest_node = (time+cost, (target,v,h))
                profit = node_sum[(time, node)] + op.functions[target](time+cost)
                if profit > node_sum[dest_node]
                    node_sum[dest_node] = profit
                    node_parent[dest_node] = (time, node) 
                end
            end
        end
    end

    #Grab the node with max profit
    score, last = findmax(node_sum)

    #Go back to reconstruct the path
    path = [last]
    parent = node_parent[last]
    while parent[1] != -1
        push!(path, parent)
        parent = node_parent[parent]
    end

    return score, reverse!(path)
end

end