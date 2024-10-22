module VsDopTvr

include("TimeFunctions.jl")
include("Helper.jl")

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
    functions::Vector{FunctionWrapper{Float64, Tuple{Float64}}}
    depots::Vector{Int64}
    tmax::Float64
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

    return OpParameters(Helper.compute_trajectories(points, graph_params), points, scores, depots, tmax)
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

function greedy_solution(op_params)
end

function variable_neighborhood_search(op_params, initial_sequence::Vector{Tuple{Int64, Int64, Int64}}, 
                                      max_iterations::Int64, verbose::Bool = false)
end


end # module VsDopTvr
