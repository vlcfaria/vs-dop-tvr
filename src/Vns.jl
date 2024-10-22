module Vns

#Operators based on
# Data Collection Planning with Non-zero Sensing Distance for a 
# Budget and Curvature Constrained Unmanned Aerial Vehicle

function waypoint_shake(sequence::Vector{Tuple{Int64,Int64,Int64}}, graph)
    #Operator randomly changes the waypoints currently used by the incumbent solution

    for i in 1:length(sequence)
        sequence[i] = (sequence[i][1], rand(1:graph.num_speeds), rand(1:graph.num_headings))
    end

    return sequence
end

function path_move(sequence::Vector{Tuple{Int64,Int64,Int64}}, graph)
    len = length(sequence)

    if len <= 1
        return sequence
    else
        start_idx = rand(2:len)
        end_idx = rand(start_idx:len)

        segment = sequence[start_idx:end_idx]
        removed = vcat(sequence[1:start_idx-1], sequence[end_idx+1:len])
        insertion_point = rand(1:len - (end_idx - start_idx + 1))

        res = vcat(removed[1:insertion_point], segment, removed[insertion_point+1:length(removed)])

        return res
    end
end

function path_exchange(sequence::Vector{Tuple{Int64,Int64,Int64}}, graph)
    len = length(sequence)
    if len < 2
        return sequence
    end

    start1 = rand(2:len)
    start2 = rand(2:len)
    #Start2 must be different
    while start2 == start1
        start2 = rand(2:len)
    end
    #Start1 should be smaller
    if start2 < start1
        start1, start2 = start2, start1
    end

    end1 = rand(start1:start2-1)
    end2 = rand(start2:len)

    seg1 = sequence[start1:end1]
    seg2 = sequence[start2:end2]

    res = vcat(sequence[1:start1-1], seg2, sequence[end1 + 1:start2-1], seg1, sequence[end2+1:len])

    return res
end

function waypoint_change(sequence::Vector{Tuple{Int64,Int64,Int64}}, graph)
    point = rand(1:length(sequence))

    sequence[point] = (sequence[point][1], rand(1:graph.num_speeds), rand(1:graph.num_headings))

    return sequence
end

function one_point_move(sequence::Vector{Tuple{Int64,Int64,Int64}}, graph)
    len = length(sequence)
    if len <= 2 #Cant apply operator since first index cannot be moved
        return sequence
    end

    idx = rand(2:len)
    val = sequence[idx]
    deleteat!(sequence, idx)

    new_idx = rand(2:len)
    #Assure index is different
    while new_idx == idx
        new_idx = rand(2:len)
    end
    insert!(sequence, new_idx, val)

    return sequence, min(new_idx, idx)
end

function one_point_exchange(sequence::Vector{Tuple{Int64,Int64,Int64}}, graph)
    len = length(sequence)
    if len <= 2 #Cant apply operator since first index cannot be moved
        return sequence
    end

    idx1 = rand(2:len)
    idx2 = rand(2:len)

    while idx1 == idx2
        idx2 = rand(2:len)
    end

    sequence[idx1], sequence[idx2] = sequence[idx2], sequence[idx1]

    return sequence, min(idx1, idx2)
end

end #end VNS