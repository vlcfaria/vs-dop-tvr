module TimeGraphHeuristic

import IterTools as itr
using DataStructures

include("AcceleratedDubins.jl")
include("Helper.jl")
include("VsDopTvr.jl")
include("TimeGraph.jl")
include("Visual.jl")

struct NodeState
    visiting_node::NTuple{3,Int64}
    visiting_timestep::Int64
    real_value::Float64
    visited_vertices::BitSet
end

#Domination relationship -> a dominates b
function Base.:>(a::NodeState, b::NodeState)
    return (a.visiting_node == b.visiting_node && a.visiting_time == b.visiting_time &&
            a.real_value >= b.real_value && issubset(a.visited_vertices, b.visited_vertices))
end

#Domination relationship -> a is dominated by b 
function Base.:<(a::NodeState, b::NodeState)
    return (a.visiting_node == b.visiting_node && a.visiting_time == b.visiting_time &&
            a.real_value < b.real_value && issubset(b.visited_vertices, a.visited_vertices))
end


function make_timestep_graphs(op, dt)
    timesteps = zeros(Int64, size(op.graph.graph))
    timesteps_to_depot = zeros(Int64, size(op.dists_to_depot))

    for i in eachindex(op.graph.graph)
        timesteps[i] = ceil(Int64, op.graph.graph[i] / dt)
    end

    for i in eachindex(op.dists_to_depot)
        timesteps_to_depot[i] = ceil(Int64, op.dists_to_depot[i] / dt)
    end

    return timesteps, timesteps_to_depot
end

#TODO - WE CAN REFINE THIS A BIT - SET TO 0 WHEN CANNOT GO BACK TO DEPOT IN TIME
function calculate_heuristic_table(op, timestep::Float64)
    steps_until_tmax = floor(Int64, op.tmax / timestep) + 1
    ans = fill(0., (length(op.functions), steps_until_tmax))

    for node in 1:length(op.functions)
        time = (steps_until_tmax - 1) * timestep
        idx = steps_until_tmax

        best = op.functions[node](time)
        while idx != 1
            idx -= 1
            time -= timestep
            best = max(best, op.functions[node](time))

            ans[node,idx] = best
        end
    end

    return ans
end

function heuristic_value(node::NodeState, table, tmax, dt)
    ans = 0.
    remaining_time = tmax - ((node.visiting_timestep - 1) * dt) #-1 to account for initial time (0, which is the timestep 1)

    for i in axes(table,1)
        if (i == node.visiting_node[1] || i in node.visited_vertices) continue end
        ans += remaining_time * table[i,node.visiting_timestep]
    end

    return ans
end

function is_pruned(f_val::Float64, cur_node::NodeState, best_node::NodeState, front_set)
    if (f_val <= best_node.real_value) return true end #1st check -> heuristic value worse than best value possible
    #Second check should be calculated in code, to prevent flooding the m_heap

    #3rd check -> check if one node in frontier is already better
    for node in front_set[cur_node.visiting_node..., cur_node.visiting_timestep]
        if issubset(node.visited_vertices, cur_node.visited_vertices) && node.real_value >= cur_node.real_value
            return true
        end
    end

    return false

end

function filter_and_add_front(cur_node::NodeState, front_set)
    idx = (cur_node.visiting_node..., cur_node.visiting_timestep)

    to_delete = NodeState[]
    for node in front_set[idx...]
        #Check if this node should be pruned
        if cur_node.real_value >= node.real_value && issubset(cur_node.visited_vertices, node.visited_vertices)
            push!(to_delete, node)
        end
    end

    for n in to_delete
        delete!(front_set, n)
    end

    push!(front_set[idx...], cur_node)
end

function heuristic_search(op, dt::Float64)
    num_nodes = length(op.functions)
    total_timesteps = floor(Int64, op.tmax / dt) + 1
    counter = 0
    starting_depot = op.depots[1] #TODO change this in the future

    heur_table = calculate_heuristic_table(op, dt)
    timestep_graph, timestep_to_depot = make_timestep_graphs(op, dt)
    m_heap = BinaryMaxHeap{Tuple{Float64, Int64, Int64, NodeState}}() #Last Int64 will be the counter, since NodeState has only partial ordering

    front_set = [Set{NodeState}() for _ in 1:num_nodes, _ in 1:op.graph.num_speeds, _ in 1:op.graph.num_headings, _ in 1:total_timesteps]
    best_node = nothing #TODO make this stable?

    for dep in op.depots
        for (v,h) in itr.product(1:op.graph.num_speeds, 1:op.graph.num_headings)
            node = NodeState((dep,v,h), 1, 0., Set([dep]))
            push!(m_heap, (node.real_value + heuristic_value(node, heur_table, op.tmax, dt), -1, counter, node))
            counter += 1
        end
    end

    while !isempty(m_heap)
        f_val, timestep, _, node = pop!(m_heap)
        timestep = -timestep
        println(f_val, node)

        if best_node !== nothing && is_pruned(f_val, node, best_node, front_set) continue end
        filter_and_add_front(node, front_set)

        if best_node === nothing || node.real_value > best_node.real_value
            best_node = node
        end

        #Define sucessors
        for target in 1:length(op.coordinates)
            #Cant visit current node and previously visited vertices
            if target == node.visiting_node[1] || target in node.visited_vertices
                continue
            end

            for (v,h) in itr.product(1:op.graph.num_speeds, 1:op.graph.num_headings)
                cost = Helper.get_dist(timestep_graph, node.visiting_node, (target,v,h))
                cost_to_depot = timestep_to_depot[target, v, h, starting_depot]

                #Valid neighbor check (respects tmax)
                if (timestep + cost + cost_to_depot > total_timesteps)
                    continue
                end
                
                profit = node.real_value + op.functions[target]((timestep+cost-1)*dt)
                new_timestep = timestep+cost
                new_visited = deepcopy(node.visited_vertices)
                push!(new_visited, target)

                dest_node = NodeState((target,v,h), new_timestep, profit, new_visited)
                heur_val = heuristic_value(dest_node, heur_table, op.tmax, dt)
                dest_f_val = heur_val + dest_node.real_value
                if is_pruned(dest_f_val, dest_node, best_node, front_set)
                    continue
                end

                push!(m_heap, (dest_f_val, -new_timestep, counter, dest_node))
                counter += 1
            end
        end


    end

    return best_node
end

end