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

function find_best_route_to_depot(op, target::Tuple{Int64,Int64,Int64})
    best_dep = 0
    best_speed = 0
    best_heading = 0

    best_time = Inf
    for (dep, speed, heading) in itr.product(1:length(op.depots), 1:op.graph.num_speeds, 1:op.graph.num_headings)
        time = op.graph.graph[target[1], dep, target[2], speed, target[3], heading]
        if time < best_time
            best_time = time
            best_dep = dep
            best_speed = speed
            best_heading = heading
        end
    end

    return (best_dep, best_speed, best_heading)
end

function calculate_seq_results(op, seq::Vector{Tuple{Int64,Int64,Int64}})
    elapsed_time = 0
    score = 0

    for i in 2:length(seq)
        prev = i-1

        #Cost from going from prev to i
        new_cost = op.graph.graph[seq[prev][1], seq[i][1], seq[prev][2], seq[i][2], seq[prev][3], seq[i][3]]

        #If cannot go back to depot in time, next target cannot be added, sequence is over
        if op.dists_to_depot[seq[i][1], seq[i][2], seq[i][3], 1] + new_cost + elapsed_time > op.tmax
            # PREV goes back to depot
            elapsed_time += op.dists_to_depot[seq[prev][1], seq[prev][2], seq[prev][3], 1]
            return score, elapsed_time, i #i is the first index not in solution
        end

        elapsed_time += new_cost
        score += op.functions[seq[i][1]](elapsed_time)
    end

    #Got to the end of vector, go back to depot anyway
    elapsed_time += op.dists_to_depot[seq[end][1], seq[end][2], seq[end][3], 1]

    return score, elapsed_time, length(seq) #all indexes are in the solution, go back
end

#Gets "actual" sequence. Includes final visit to depot, and stops when visiting final visit.
function get_actual_sequence(op, seq::Vector{Tuple{Int64, Int64, Int64}})
    elapsed_time = 0
    full_seq = Vector{Tuple{Int64, Int64, Int64}}()
    push!(full_seq, seq[1]) #push depot

    for i in 2:length(seq)
        prev = i - 1
        new_cost = op.graph.graph[seq[prev][1], seq[i][1], seq[prev][2], seq[i][2], seq[prev][3], seq[i][3]]

        #Check if can add point and still return to depot
        if op.dists_to_depot[seq[i][1], seq[i][2], seq[i][3], 1] + new_cost + elapsed_time > op.tmax
            push!(full_seq, find_best_route_to_depot(op, seq[prev]))

            return full_seq
        end

        #Can continue
        elapsed_time += new_cost
        push!(full_seq, seq[i])
    end

    push!(full_seq, find_best_route_to_depot(op, seq[end]))
end

function retrieve_path(op_params, configurations::Vector{Tuple{Int64, Int64, Int64}})
    #(dubinspath, starting_speed, ending_speed)
    full_path::Vector{Tuple{AcceleratedDubins.DubinsPathR2, Float64, Float64}} = []
    configurations = get_actual_sequence(op_params, configurations)

    locations = op_params.coordinates
    headings = op_params.graph.headings
    speeds = op_params.graph.speeds
    radii = op_params.graph.radii

    vehicle_params = op_params.graph.vehicle_params
    params = [vehicle_params.v_min, vehicle_params.v_max, vehicle_params.a_max, -vehicle_params.a_min]

    for i in 2:length(configurations)
        prev = i - 1

        n_i, v_i, h_i = configurations[prev]
        n_f, v_f, h_f = configurations[i]

        starting::Vector{Float64} = [locations[n_i][1], locations[n_i][2], headings[h_i]]
        ending::Vector{Float64} = [locations[n_f][1], locations[n_f][2], headings[h_f]]

        path, _, _ = AcceleratedDubins.fastest_path(starting, ending, radii, params, [speeds[v_i], speeds[v_f]])
        push!(full_path, (path, speeds[configurations[prev][2]], speeds[configurations[i][2]]))
    end

    return full_path
end

function compute_trajectories(locations::Vector{Tuple{Float64, Float64}}, graph_params)
    #Build a 6-dimensional graph (starting node, ending node, starting speed, ending speed, starting heading angle, ending heading angle)
    num_locations = length(locations)
    num_speeds = graph_params.num_speeds
    num_headings = graph_params.num_headings
    
    graph = graph_params.graph

    speeds = graph_params.speeds
    headings = graph_params.headings
    radii = graph_params.radii

    # v_min v_max a_max -a_min
    vehicle_params = graph_params.vehicle_params
    params = [vehicle_params.v_min, vehicle_params.v_max, vehicle_params.a_max, -vehicle_params.a_min]
    
    @Threads.threads for node_i in 1:num_locations
        for node_f in 1:num_locations
            for (v_i, v_f) in itr.product(1:num_speeds, 1:num_speeds)
                for (h_i, h_f) in itr.product(1:num_headings, 1:num_headings)
                    start::Vector{Float64} = [locations[node_i][1], locations[node_i][2], headings[h_i]]
                    stop::Vector{Float64} = [locations[node_f][1], locations[node_f][2], headings[h_f]]
                    path, time, _ = AcceleratedDubins.fastest_path(start, stop, radii, params, [speeds[v_i], speeds[v_f]])
                    #Set to edge, note that we can't take advantage of simmetry because acceleration max/min is not necessarily the same
                    graph[node_i,node_f,v_i,v_f,h_i,h_f] = path === nothing ? Inf : time
                end
            end
        end
    end

    return graph_params
end

function compute_fastests_paths_to_depot(graph, depots)
    num_targets = size(graph.graph)[1]
    num_speeds = graph.num_speeds
    num_headings = graph.num_headings

    fastest_paths = Array{Float64, 4}(undef, num_targets, num_speeds, num_headings, length(depots))

    for dep in 1:length(depots)
        for (target, speed, heading) in itr.product(1:num_targets, 1:num_speeds, 1:num_headings)
            # Check for fastest path in graph
            best = Inf
            for (dep_speed, dep_heading) in itr.product(1:num_speeds, 1:num_headings)
                time = graph.graph[target, depots[dep], speed, dep_speed, heading, dep_heading]
                if time < best
                    best = time
                end
            end

            #Fill in
            fastest_paths[target, speed, heading, dep] = best
        end
    end

    return fastest_paths
end

end