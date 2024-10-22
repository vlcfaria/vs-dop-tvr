module Helper

function read_optvr_file(filename::String)
    coordinates::Array{Tuple{Float64, Float64}, 1} = []
    values::Vector{Float64} = []
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
            elseif line == "NODE_SCORE_SECTION"
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
                score = parse(Float64, tokens[2])
                push!(values, score)
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

end