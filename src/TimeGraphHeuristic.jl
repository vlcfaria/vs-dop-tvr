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
    visited_vertices::Set{Int64}
end

#Domination relationship
function Base.:>(a::NodeState, b::NodeState)
    return (a.visiting_node == b.visiting_node && a.visiting_time == b.visiting_time &&
            a.real_value >= b.real_value && issubset(a.visited_vertices, b.visited_vertices))
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
    remaining_time = tmax - ((node.visiting_timestep - 1) * dt) #-1 to account for initial time (0)

    for i in 1:axes(table,1)
        if (i in node.visited_vertices) continue end
        ans += remaining_time * table[i,node.visiting_timestep]
    end

    return ans
end

function is_pruned(f_val::Float64, cur_node::NodeState, best_node::NodeState, front_set)
    if (f_val <= best_node.real_value) return true end #1st check -> heuristic value worse than best value possible
    #Second check should be calculated in code, to prevent flooding the m_heap

    #3rd check -> check if one node in frontier is already better
    for node in front_set[cur_node.visiting_node, cur_node.visiting_timestep]
        if issubset(node.visited_vertices, cur_node.visited_vertices) && node.real_value >= cur_node.real_value
            return true
        end
    end

    return false

end

function filter_and_add_front(cur_node::NodeState, front_set)
    for node in front_set[node.visiting_node, node.visiting_timestep]
        #Check if this node should be pruned
        if cur_node.real_value >= node.real_value && issubset(cur_node.visited_vertices, node.visited_vertices)
            delete!(front_set, node)
        end
    end

    push!(front_set, node)
end

function heuristic_search(op, dt::Float64)
    num_nodes = length(op.functions)
    total_timesteps = floor(Int64, op.tmax / dt) + 1

    heur_table = calculate_heuristic_table(op, dt)
    m_heap = BinaryMaxHeap{Tuple{Float64, Int64, NodeState}}()
    #Also calculate graph in number of timesteps between two nodes

    front_set = [Set{NodeState}() for i in 1:num_nodes for j in 1:total_timesteps]
    best_node = nothing #TODO make this stable?

    for dep in op.depots
        for (v,h) in itr.product(1:op.graph.num_speeds, 1:op.graph.num_headings)
            node = NodeState((dep,v,h), 1, 0., Set([dep]))
            push!(m_heap, (node.real_value + heuristic_value(node, heur_table, op.tmax, dt), 1, node))
        end
    end

    while !isempty(m_heap)
        f_val, timestep, node = pop!(m_heap)

        if is_pruned(f_val, node, best_node, front_set) continue end
        filter_and_add_front(node, front_set)


    end

end

end