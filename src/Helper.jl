module Helper

include("TimeFunctions.jl")

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

end