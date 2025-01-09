module Vns

using Random

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
        
        #Change all waypoints in the segment
        segment = [(p, rand(1:graph.num_speeds), rand(1:graph.num_headings)) for (p,_,_) in segment]

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

    #Change waypoints of both segments
    seg1 = [(p, rand(1:graph.num_speeds), rand(1:graph.num_headings)) for (p,_,_) in seg1]
    seg2 = [(p, rand(1:graph.num_speeds), rand(1:graph.num_headings)) for (p,_,_) in seg2]

    res = vcat(sequence[1:start1-1], seg2, sequence[end1 + 1:start2-1], seg1, sequence[end2+1:len])

    return res
end

function waypoint_change(sequence::Vector{Tuple{Int64,Int64,Int64}}, graph, limit_idx)
    point = rand(1:limit_idx) #Only change up to limit index, changing others wont do nothing new

    sequence[point] = (sequence[point][1], rand(1:graph.num_speeds), rand(1:graph.num_headings))

    return sequence, point
end

function one_point_move(sequence::Vector{Tuple{Int64,Int64,Int64}}, graph, limit_idx)
    len = length(sequence)
    if len <= 2 #Cant apply operator since first index cannot be moved
        return sequence
    end

    idx = rand(1:len)
    val = sequence[idx]
    deleteat!(sequence, idx)

    if idx == 1 #If changing first point (depot), keep it there
        new_idx = 1
    elseif idx <= limit_idx #Any chosen index will modify the sequence
        new_idx = rand(2:len)
        #Assure index is different
        while new_idx == idx
            new_idx = rand(2:len)
        end
    else #Chosen index is out of the limit index, new index must be up to limit
        new_idx = rand(2:limit_idx)
    end

    #Randomly change waypoint too
    val = (val[1], rand(1:graph.num_speeds), rand(1:graph.num_headings))
    insert!(sequence, new_idx, val)

    return sequence, min(new_idx, idx)
end

function one_point_exchange(sequence::Vector{Tuple{Int64,Int64,Int64}}, graph, limit_idx)
    len = length(sequence)
    if len <= 2 #Cant apply operator since first index cannot be moved
        return sequence
    end

    idx1 = rand(2:len)

    if idx1 <= limit_idx
        idx2 = rand(2:len) #Choosing anything will modify the sequence
        while idx1 == idx2 #While also assuring it is different
            idx2 = rand(2:len)
        end
    else
        #Need to choose a value within limit index
        idx2 = rand(2:limit_idx)
    end

    #Swap and assign random waypoints
    sequence[idx1], sequence[idx2] = sequence[idx2], sequence[idx1]

    sequence[idx1] = (sequence[idx1][1], rand(1:graph.num_speeds), rand(1:graph.num_headings))
    sequence[idx2] = (sequence[idx2][1], rand(1:graph.num_speeds), rand(1:graph.num_headings))

    return sequence, min(idx1, idx2)
end

function shake(sequence::Vector{Tuple{Int64,Int64,Int64}}, graph, l::Int64)
    if l == 1
        return waypoint_shake(sequence, graph)
    elseif l == 2
        return path_move(sequence, graph)
    else
        return path_exchange(sequence, graph)
    end
end


#limit_idx is the first index of the solution that is not visited, or the last index, in case all points is visited
function search(sequence::Vector{Tuple{Int64,Int64,Int64}}, graph, l::Int64, limit_idx::Int64)
    if l == 1
        return waypoint_change(sequence, graph, limit_idx)
    elseif l == 2
        return one_point_move(sequence, graph, limit_idx)
    else
        return one_point_exchange(sequence, graph, limit_idx)
    end
end

end #end VNS