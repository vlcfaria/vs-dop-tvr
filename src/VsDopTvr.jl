module VsDopTvr

import IterTools as itr
using Random

include("TimeFunctions.jl")
include("AcceleratedDubins.jl")
include("Helper.jl")
include("Visual.jl")
include("Vns.jl")
include("LocalSearch.jl")

using FunctionWrappers
import FunctionWrappers: FunctionWrapper


struct VehicleParameters
    v_min::Float64
    v_max::Float64
    a_min::Float64
    a_max::Float64
    r_min::Float64
    r_max::Float64
end

struct DOPGraph
    graph::Array{Float64, 6}
    vehicle_params::VehicleParameters
    num_speeds::Int64
    num_headings::Int64
    radii_samples::Int64
    speeds::Vector{Float64}
    headings::Vector{Float64}
    radii::Vector{Float64}

    function DOPGraph(num_locations::Int64, v_params)
        num_speeds = 4
        num_headings = 8
        radii_samples = 8

        return DOPGraph(num_speeds, num_headings, radii_samples, num_locations, v_params)
    end

    function DOPGraph(speeds::Vector{Float64}, headings::Vector{Float64}, 
                      radii_samples::Int64, num_locations::Int64, v_params)
        num_speeds = length(speeds)
        num_headings = length(headings)
        radii = AcceleratedDubins.radii_samples_exp(v_params.r_min, v_params.r_max, radii_samples)
        graph = fill(-1, (num_locations, num_locations, num_speeds, num_speeds, num_headings, num_headings))

        return new(graph, v_params, num_speeds, num_headings, radii_samples, deepcopy(speeds), deepcopy(headings), radii)
    end

    function DOPGraph(num_speeds::Int64, num_headings::Int64, 
                      radii_samples::Int64, num_locations::Int64, v_params)
        
        # Max discrete value at curve
        max_speed_discrete = min(v_params.v_max, AcceleratedDubins.speed_by_radius(v_params.r_max))
        speeds = collect(range(v_params.v_min, max_speed_discrete, num_speeds))
        headings = collect(range(0, 2 * pi, num_headings))
        radii = AcceleratedDubins.radii_samples_exp(v_params.r_min, v_params.r_max, radii_samples)
        graph = fill(-1, (num_locations, num_locations, num_speeds, num_speeds, num_headings, num_headings))

        return new(graph, v_params, num_speeds, num_headings, 
                   radii_samples, speeds, headings, radii)
    end
end

mutable struct OpParameters
    graph::DOPGraph
    coordinates::Array{Tuple{Float64, Float64}, 1}
    functions::Vector{TimeFunctions.TimeFunc}
    depots::Vector{Int64}
    tmax::Float64
    dists_to_depot::Array{Float64, 4} #TODO this can be a 3-dimensional array, calculate to fastest depot in case of more depots
    #(n,v,h,depot)
end

mutable struct Results
    scores::Vector{Float64}
    travel_times::Vector{Float64}
    timestamps::Vector{Float64}
end 

function get_cessna172_params()
    return VehicleParameters(30., 67., -3., 2., 65.7, 264.2)
end

function build_graph(instance_path::String, vehicle_params, graph_params = nothing)
    points, scores, depots, tmax = Helper.read_optvr_file(instance_path)
    
    if graph_params === nothing
        graph_params = DOPGraph(length(points), vehicle_params)
    end

    graph_params = Helper.compute_trajectories(points, graph_params)
    return OpParameters(graph_params, points, scores, depots, tmax, Helper.compute_fastests_paths_to_depot(graph_params, depots))
end

function vs_dop_tvr(op_params, max_iterations::Int64 = 2000, verbose::Bool = false)
    initial_seq, ini_time, ini_score = greedy_solution(op_params)
    
    if verbose
        println("Initial Score: ", ini_score, ". Initial time: ", ini_time)
    end

    vns_seq, time, score = variable_neighborhood_search(op_params, initial_seq, max_iterations, verbose)

    config = Helper.shortest_configuration_by_sequence(op_params, vns_seq)
    path = Helper.retrieve_path(op_params, config)

    return path, config, score, time
end

function greedy_solution(op)
    to_add = Set{Int64}(1:length(op.coordinates))

    #Remove all depots as candidates, since they dont give score
    for d in op.depots
        delete!(to_add, d)
    end

    num_speeds = op.graph.num_speeds
    num_headings = op.graph.num_headings

    #Select the best first choice for depot, including speed & heading
    starting_pos = (-1,-1,-1)
    second_pos = (-1,-1,-1)
    best_ratio = 0

    for dep in op.depots #For all depots
        for cand in to_add # Check all candidates
            for (v_d, h_d, v_s, h_s) in itr.product(1:num_speeds, 1:num_headings, 1:num_speeds, 1:num_headings)
                time_inc = op.graph.graph[dep, cand, v_d, v_s, h_d, h_s]
                score_inc = op.functions[cand](time_inc)
                ratio = score_inc / time_inc

                #If move is allowed and ratio is better
                if ratio > best_ratio && time_inc + op.dists_to_depot[cand, v_s, h_s, dep] <= op.tmax
                    best_ratio = ratio
                    starting_pos = (dep, v_d, h_d)
                    second_pos = (cand, v_s, h_s)
                end
            end
        end
    end

    sequence::Vector{Tuple{Int64,Int64,Int64}} = [starting_pos, second_pos]
    delete!(to_add, second_pos[1])

    seq_time = op.graph.graph[starting_pos[1], second_pos[1], starting_pos[2], second_pos[2], starting_pos[3], second_pos[3]]
    seq_score = op.functions[second_pos[1]](seq_time)

    can_add = true
    while can_add
        can_add = false

        prev = sequence[end]
        best_score_increment = 0
        best_ratio = 0
        best_time_increment = 0
        best_move = (-1,-1,-1)

        # Find next best move
        for candidate in to_add
            # Check every edge and its legality (can go back to depot)
            for (v,h) in itr.product(1:num_speeds, 1:num_headings)
                time_increment = op.graph.graph[prev[1], candidate, prev[2], v, prev[3], h]
                score_increment = op.functions[candidate](seq_time + time_increment)
                ratio = score_increment / time_increment
    
                if ratio > best_ratio && op.dists_to_depot[candidate, v, h, 1] + time_increment + seq_time <= op.tmax
                    best_move = (candidate, v, h)
                    best_score_increment = score_increment
                    best_time_increment = time_increment
                    best_ratio = ratio
                    can_add = true
                end
            end
        end
        
        if can_add #append best move
            delete!(to_add, best_move[1])
            push!(sequence, best_move)
            seq_time += best_time_increment
            seq_score += best_score_increment
        end
    end

    last = sequence[end]
    seq_time += op.dists_to_depot[last[1], last[2], last[3]]

    remaining = shuffle(collect(to_add))
    sequence = vcat(sequence, [(n, rand(1:num_speeds), rand(1:num_headings)) for n in remaining])

    return sequence, seq_score, seq_time
end

function random_solution(op)
    starting_depot = rand(op.depots)

    #Get nodes and remove chosen depot
    nodes = collect(1:length(op.coordinates))
    deleteat!(nodes, starting_depot)
    shuffle!(nodes)
    insert!(nodes, 1, starting_depot)

    speeds, headings = op.graph.num_speeds, op.graph.num_headings
    seq = [(n, rand(1:speeds), rand(1:headings)) for n in nodes] #Create random solution

    score, time, _= Helper.calculate_seq_results(op, seq)
    return seq, score, time
end

function variable_neighborhood_search(op, initial_sequence::Vector{Tuple{Int64, Int64, Int64}}, 
                                      max_iterations::Int64, verbose::Bool = false)

    best_score, best_time, _ = Helper.calculate_seq_results(op, initial_sequence)
    best_sequence = initial_sequence
    len = length(initial_sequence)
    if (best_time > op.tmax) 
        throw("Invalid initial solution given")
    end

    for i in 1:max_iterations
        if verbose
            println(("OUTER", i, best_score, best_time))
        end

        if i % 250 == 0
            println((i, best_score, best_time, Helper.get_actual_sequence(op, best_sequence)))
        end

        l = 1
        while l <= 3
            #Shake
            local_sequence = Vns.shake(deepcopy(best_sequence), op.graph, l)
            #TODO useless line below, remove
            local_score, local_time, local_limit_idx = Helper.calculate_seq_results(op, local_sequence)

            if verbose
                println(("Start", l, local_score, local_time))
            end
            
            #Search
            local_seq, local_score, local_time = LocalSearch.search(deepcopy(local_sequence), op, l)

            if verbose
                println(("Final", local_score, local_time))
            end
            
            #Higher score found through local search OR equal score with lower time
            if (local_time <= op.tmax && local_score > best_score) || (local_score == best_score && local_time < best_time)
                best_time = local_time
                best_sequence = local_seq
                best_score = local_score
                l = 1
            else
                l += 1
            end
        end
    end

    return best_sequence, best_time, best_score
end

function general_vns(op, seq, max_iter, l_max=3)
    score, time, _ = Helper.calculate_seq_results(op,seq)

    #General VNS algorithm
    for _ in 1:max_iter
        l = 1
        while l <= l_max
            #Shake
            l_seq = Vns.shake(deepcopy(seq), op.graph, l)

            #Call variable neighborhood descent (VND)
            l_seq, l_score, l_time = LocalSearch.VND(l_seq, op, l_max)

            #Change neighborhood up or down
            if l_score > score || l_score == score && l_time < time
                seq, score, time = l_seq, l_score, l_time
                l = 1
            else
                l += 1
            end
        end
    end

    return seq, score
end

end # module VsDopTvr
