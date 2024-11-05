module Helper

include("TimeFunctions.jl")
import IterTools as itr
using ..AcceleratedDubins

function read_optvr_file(filename::String)
    coordinates::Array{Tuple{Float64, Float64}, 1} = []
    values::Vector{TimeFunctions.TimeFunc} = []
    depots::Vector{Int64} = []
    tmax = -1.

    open(filename) do file
        in_coordinates = false
        in_scores = false
        in_depots = false
        
        for line in eachline(file)
            tokens = split(line)
            if line == "NODE_COORD_SECTION"
                in_coordinates = true
            elseif line == "NODE_FUNC_SECTION"
                in_coordinates = false
                in_scores = true
            elseif line == "DEPOT_SECTION"
                in_scores = false
                in_depots = true
            elseif tokens[1] == "COST_LIMIT"
                tmax = parse(Float64, tokens[3])
            elseif line == "EOF"
                break
            elseif in_coordinates
                x_coord = parse(Float64, tokens[2])
                y_coord = parse(Float64, tokens[3])
                push!(coordinates, (x_coord, y_coord))
            elseif in_scores
                #Parse rest
                args = [parse(Float64, x) for x in tokens[3:end]]
                #TODO make this a map
                if (tokens[2] == "sine_wave")
                    func = TimeFunctions.sine_wave(args...)
                else
                    throw(ArgumentError("Function not recognized: " * tokens[2]))
                end
                push!(values, func)
            elseif in_depots
                depot = parse(Int64, tokens[1])
                if depot == -1
                    in_depots = false
                else
                    push!(depots, depot)
                end                  
            end
        end
    end

    return coordinates, values, depots, tmax
end

function is_legal_move(op, seq::Vector{Tuple{Int64,Int64,Int64}}, depot, time)
    #Check if last can return
    last = seq[end]
    for (v,h) in itr.product(1:op.graph.num_speeds, 1:op.graph.num_headings)
        if time + op.graph.graph[last[1], op.depots[1], last[2], v, last[3], h] <= op.tmax
            return true
        end
    end
    return false
end

function calculate_seq_score(op, seq::Vector{Tuple{Int64,Int64,Int64}})
    elapsed_time = 0
    score = 0

    for i in 2:length(seq)
        if seq[i][1] == op.depots[1]
            break
        end
        prev = i-1

        elapsed_time += op.graph.graph[seq[prev][1], seq[i][1], seq[prev][2], seq[i][2], seq[prev][3], seq[i][3]]
        score += op.functions[seq[i][1]](elapsed_time)
    end

    return score, elapsed_time
end

function retrieve_path(op_params, configurations::Vector{Tuple{Int64, Int64, Int64}})
    #(dubinspath, starting_speed, ending_speed)
    full_path::Vector{Tuple{AcceleratedDubins.DubinsPathR2, Float64, Float64}} = []

    locations = op_params.coordinates
    headings = op_params.graph.headings
    speeds = op_params.graph.speeds
    radii = op_params.graph.radii

    vehicle_params = op_params.graph.vehicle_params
    params = [vehicle_params.v_min, vehicle_params.v_max, vehicle_params.a_max, -vehicle_params.a_min]

    for i in 1:length(configurations)-1
        if i != 1 && configurations[i][1] == op_params.depots[1]
            break
        end
        
        next = i + 1

        n_i, v_i, h_i = configurations[i]
        n_f, v_f, h_f = configurations[next]

        starting::Vector{Float64} = [locations[n_i][1], locations[n_i][2], headings[h_i]]
        ending::Vector{Float64} = [locations[n_f][1], locations[n_f][2], headings[h_f]]

        path, time, _ = AcceleratedDubins.fastest_path(starting, ending, radii, params, [speeds[v_i], speeds[v_f]])

        push!(full_path, (path, speeds[configurations[i][2]], speeds[configurations[next][2]]))
    end

    return full_path
end

end