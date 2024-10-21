## ADAPTED FROM https://github.com/comrob/AcceleratedDubins.jl

"""BSD 2-Clause License

Copyright (c) 2020, Computational Robotics Laboratory 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."""

module AcceleratedDubins

###############################################################################
############################ TYPE DEFINITIONS #################################
###############################################################################

@enum DubinsPathTypeR2 LSL LSR RSL RSR
@enum SegmentType L_SEG S_SEG R_SEG

DIRDATA = Dict{Int,Vector{SegmentType}}(
    Int(LSL) => [L_SEG, S_SEG, L_SEG],
    Int(LSR) => [L_SEG, S_SEG, R_SEG],
    Int(RSL) => [R_SEG, S_SEG, L_SEG],
    Int(RSR) => [R_SEG, S_SEG, R_SEG],
)

const PATH_OK = 0       # no error, path found
const PATH_ERROR = 1    # path not possible

"""
The data structure that holds the full dubins path with multiple radii.

# Data fields
- `origin::Vector{Float64}`: the initial configuration [x, y, a], angle in radians
- `lengths::Vector{Float64}`: lengths of each segment
- `r::Vector{Float64}`: turn radii of each segment, where radius of straight path is Inf
- `type::DubinsPathTypeR2`: the Dubins path type given by the @enum DubinsPathTypeR2
"""
mutable struct DubinsPathR2
   origin::Vector{Float64}
   lengths::Vector{Float64}
   r::Vector{Float64}
   type::DubinsPathTypeR2
end

DubinsPathR2() = DubinsPathR2(zeros(3), zeros(3), [1, Inf, 1], LSL)

###############################################################################
########################### HELPER FUNCTIONS ##################################
###############################################################################
"""
    radii_samples_exp(r_min::Float64, r_max::Float64, num_samples::Int64, precision::Int64=100)
Calculate radii samples with exponential spacing.
...
# Arguments
- `r_min::Float64`: minimum turning radius
- `r_max::Float64`: maximum turning radius
- `num_samples::Int64`: number of samples
- `precision::Int64`: for rounding the values. Default: Inf (no rounding)
...
return - samples::Vector{Float64}

See also: [`radii_samples_lin`](@ref)
"""
function radii_samples_exp(r_min::Float64, r_max::Float64, num_samples::Int64, precision::Int64=100)
    if precision == 100
        return [exp(i) for i in log(r_min):(log(r_max) - log(r_min))/(num_samples-1):log(r_max)]
    end
    return [round(exp(i), digits=precision) for i in log(r_min):(log(r_max) - log(r_min))/(num_samples-1):log(r_max)]
end

" Compute DubinsR2 path length. Argument: `path::DubinsPathR2, return: `len::Float64`"
function path_len(path)
    return (path === nothing) ? Inf : sum(path.lengths)
end

"""
    is_path_faster(path_lengths, radii::Vector{Float64}, best_time::Float64, params, speeds)

Comparing function to use between two maneuvres.
...
# Arguments
- `path_lengths::Vector{Float64}`: vector of segment lengths
- `radii::Vector{Float64}`: vector of segment radii
- `best_time::Float64`: time to compare with
- `params::Vector{Float64}`: Vector of speed parameters: [velocity in minimal radius, maximal vehicle velocity, acceleration, deceleration]
- `speeds::Vector{Float64}` : Vector of speed values: [starting velocity, ending velocity]
...
return - ans::Boolean, time::Float64
"""
function is_path_faster(path_lengths, radii::Vector{Float64}, best_time::Float64, params, speeds)
    if path_lengths === nothing
        return false, nothing
    end
    path = DubinsPathR2()
    path.type = LSL # for time computing, only posible maneuvres currently are CSC
    path.lengths = path_lengths
    path.r = radii
    time = path_time(path, params, speeds)
    if time < best_time
        return true, time
    end
    return false, time
end

"""
    speed_by_radius(radius::Number)
Caculation using plane turn

# Arguments
- `radius::Number`: Radius to get speed for

# Return
- `v::Float64`: Maximum speed for input radius
"""
function speed_by_radius(radius::Number)
    v = sqrt(radius*tan(pi/3)*9.80665) #g
    return v
end

"Opposite function"
function radius_by_speed(speed::Number)
    r = speed^2 / (tan(pi/3)*9.80665)
    return r
end

"""
    path_time(path::DubinsPathR2, params::Vector{Number}, speeds::Vector{Float64})
Get path time of single path.

# Arguments
- `path::DubinsPathR2`: Path to compute time from
- `params::Vector{Float64}`: Vector of speed parameters: [velocity in minimal radius, maximal vehicle velocity, acceleration, deceleration]
- `speeds::Vector{Float64}` : Vector of speed values: [starting velocity, ending velocity]

# Return
- `time::Float64`: Path time
"""
function path_time(path::DubinsPathR2, params::Vector{Float64}, speeds::Vector{Float64})
    if path === nothing
        return Inf
    end
    t, _ = speed_profile(path, params[1], params[2], params[3], params[4], speeds[1], speeds[2])
    return t[end]
end

"""
    point_by_angle(origin::Vector{Float64}, ang::Number, len::Number)
Get point that is in certain distance by certain angle from origin point

# Arguments
- `origin::Vector{Float64}`: starting angle in radians from point [0, 0]
- `ang::Number`: starting angle in radians from point [0, 0]
- `len::Number`: starting angle in radians from point [0, 0]

# Return
- `[x2,y2]::Vector{Float64}`: desired point
"""
function point_by_angle(origin::Vector{Float64}, ang::Number, len::Number)
    x, y = origin[1], origin[2]
    y2 = y + len * sin(ang)
    x2 = x + len * cos(ang)
    return [x2, y2]
end

"""
    rotate_points(start::Vector{Float64}, stop::Vector{Float64})
Rotate given input points to format start -> [0, 0, alfa], stop -> [dist, 0, beta]

# Arguments
- `start::Vector{Float64}`: start point of the vehicle, [x1, y1, theta1]
- `stop::Vector{Float64}`: final point of the vehicle, [x2, y2, theta2]

# Return
- `alfa::Float64`: new starting angle
- `beta::Float64`: new ending angle
- `dist::Float64`: distance between start and stopex
"""
function rotate_points(start::Vector{Float64}, stop::Vector{Float64})
    diff = stop-start
    global_rotation = atan(diff[2], diff[1])
    alfa = start[3] - global_rotation
    beta = stop[3] - global_rotation
    dist = sqrt(diff[1]^2 + diff[2]^2)
    return alfa, beta, dist
end

"""
    function possible_tuples(radii::Vector{Float64})
Get all possible combinations of input radii

# Arguments
- `radii::Vector{Float64}`: array of possible radii

# Return
- `tup::Vector{Float64}`: pairs of radii
"""
function possible_tuples(radii::Vector{Float64})
    tup::Vector{Vector{Float64}} = []
    for i in 1  :length(radii)
        for j in i:length(radii)
            push!(tup, [radii[i], radii[j]])
            if i != j
                push!(tup, [radii[j], radii[i]])
            end
        end
    end
    return tup
end

###############################################################################
######################### MAIN PATH FUNCTIONS #################################
###############################################################################
"""
    fastest_path(start::Vector{Float64}, stop::Vector{Float64}, radii::Vector{Float64}, params::Array{Float64}, speeds::Array{Float64})
based on the input radii, tests each combination and select best DubinsR2 path.

# Arguments
- `start::Vector{Float64}`: starting angle in radians from point [0, 0]
- `stop::Vector{Float64}`: starting angle in radians from point [0, 0]
- `radii::Vector{Float64}`: radius of starting arc
- `params::Array{Float64}`: Vector of speed parameters: [velocity in minimal radius, maximal vehicle velocity, acceleration, deceleration]
- `speeds::Vector{Float64}` : Vector of speed values: [starting velocity, ending velocity]

# Return
- `path::DubinsPathR2`: fastest DubinsR2 path
- `err::const`: error code of the path calculation
"""
function fastest_path(start::Vector{Float64}, stop::Vector{Float64}, radii::Vector{Float64}, params::Array{Float64}, speeds::Array{Float64})
    alfa, beta, dist = rotate_points(start, stop)

    # get all possible combinations + best possible radius by the speeds
    radii_pairs = push!(possible_tuples(radii), [radius_by_speed(speeds[1]), radius_by_speed(speeds[2])])

    best_time = Inf
    best_len = Inf
    best_path = nothing
    ret_err = PATH_ERROR
    for pair in radii_pairs
        path, err = best_maneuver(alfa, beta, dist, pair[1], pair[2], params, speeds)
        if path === nothing
            continue
        end
        time = path_time(path, params, speeds)
        if time < best_time || (isapprox(time, best_time, atol=1e-8) && path_len(path) < best_len)
            best_time = time
            best_len = path_len(path)
            best_path = path
            ret_err = err
        end
    end

    if best_path !== nothing
        best_path.origin = start
    end
    return best_path, best_time, ret_err
end


"""
    best_maneuver(alfa::Float64, beta::Float64, dist::Float64, r1::Float64, r2::Float64, params::Array{Float64}, speeds::Array{Float64}))
Function to return best DubinsR2 path according to is_path_better argument

# Arguments
- `alfa::Vector{Float64}`: starting angle in radians from point [0, 0]
- `beta::Float64`: ending angle in radians from point [dist, 0]
- `dist::Float64`: distance between starting and ending point
- `r1::Float64`: radius of starting arc
- `r2::Float64`: radius of ending arc
- `params::Array{Float64}`: Vector of speed parameters: [velocity in minimal radius, maximal vehicle velocity, acceleration, deceleration]
- `speeds::Vector{Float64}` : Vector of speed values: [starting velocity, ending velocity]

# Return
- `path::DubinsPathR2`: path 
- `err`: result of operation
"""
function best_maneuver(alfa::Float64, beta::Float64, dist::Float64, r1::Float64, r2::Float64, params::Array{Float64}, speeds::Array{Float64})
    path = DubinsPathR2()
    path.origin = [0, 0, alfa]
    best = Inf

    "shortcut function to set all the path variables"
    function pset(flengths, fbest, fr1, fr2, fr3, ftype)
        path.lengths = flengths
        best = fbest
        path.r[1] = fr1
        path.r[2] = fr2
        path.r[3] = fr3
        path.type = DubinsPathTypeR2(ftype)
    end

    _, path_lengths = RSR_path(alfa, beta, dist, r1, r2)
    ans, calc = is_path_faster(path_lengths, [r1, Inf, r2], best, params, speeds)
    if ans
        pset(path_lengths, calc, r1, Inf, r2, 3)
    end

    _, path_lengths = LSL_path(alfa, beta, dist, r1, r2)
    ans, calc = is_path_faster(path_lengths, [r1, Inf, r2], best, params, speeds)
    if ans
        pset(path_lengths, calc, r1, Inf, r2, 0)
    end

    _, path_lengths = RSL_path(alfa, beta, dist, r1, r2)
    ans, calc = is_path_faster(path_lengths, [r1, Inf, r2], best, params, speeds)
    if ans
        pset(path_lengths, calc, r1, Inf, r2, 2)
    end

    _, path_lengths = LSR_path(alfa, beta, dist, r1, r2)
    ans, calc = is_path_faster(path_lengths, [r1, Inf, r2], best, params, speeds)
    if ans
        pset(path_lengths, calc, r1, Inf, r2, 1)
    end

    if best == Inf
        return nothing, PATH_ERROR
    end
    return path, PATH_OK
end

###############################################################################
#################### PATH SAMPLING FOR PLOTTING ###############################
###############################################################################
"""
    sample_path(path::DubinsPathR2, resolution::Float64 = 0.01)
Sample single path for plotting

# Arguments
- `path::DubinsPathR2`: path to sample
- `resolution::Float64`: how smoothly to sample the arc, smaller -> finer

# Return
- `x::Vector{Float64}`: x points of the path to plot
- `y::Vector{Float64}`: y points of the path to plot
"""
function sample_path(path::DubinsPathR2, resolution::Float64 = 0.01)
    if path === nothing
        return Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
    else
        confx::Vector{Vector{Float64}} = [[], [], []]
        confy::Vector{Vector{Float64}} = [[], [], []]

        # Sample first arc
        for l in 0.:resolution:path.lengths[1]
            x, y, _ = get_configuration(path, l)
            push!(confx[1], x)
            push!(confy[1], y)
        end

        # Add last point of first arch
        x, y, theta = get_configuration(path, path.lengths[1])
        push!(confx[2], x)
        push!(confy[2], y)

        # Sample line segment
        x, y = point_by_angle([x, y], theta, path.lengths[2])       
        push!(confx[2], x)
        push!(confy[2], y)

        # Sample third segment
        s3_start = path.lengths[1] + path.lengths[2]
        s3_end = sum(path.lengths)
        for l in s3_start:resolution:s3_end
            x, y, _ = get_configuration(path, l)
            push!(confx[3], x)
            push!(confy[3], y)
        end

        # Add last point of third arch
        x, y, theta = get_configuration(path, s3_end)
        push!(confx[3], x)
        push!(confy[3], y)

        return confx, confy
    end
end


###############################################################################
###################### PATH TYPES COMPUTATION #################################
###############################################################################

function RSR_path(alfa::Float64, beta::Float64, dist::Number, radius1::Number, radius2::Number)
    S1 = (radius1*sin(alfa), -radius1*cos(alfa))
    S2 = (dist+radius2*sin(beta), -radius2*cos(beta))

    dx = S2[1] - S1[1]
    dy = S2[2] - S1[2]
    d = sqrt(dx^2 + dy^2)
    if (d^2 - (radius1 - radius2)^2) < 0
        return PATH_ERROR, nothing
    end
    l = sqrt(d^2 - (radius1 - radius2)^2)

    gamma = atan(dy/dx)
    xi = acos((radius1-radius2)/d)

    mu = mod2pi(gamma+xi)
    if (S2[1] < S1[1])
        mu += pi
    end

    out = Vector{Float64}(undef,3)
    out[1] = mod2pi(alfa + pi/2 - mu) * radius1
    out[2] = l
    out[3] = mod2pi(2*pi - (beta + pi/2 - mu)) * radius2

    return PATH_OK, out
end

function LSL_path(alfa, beta, dist, radius1, radius2)
    S1 = (-radius1*sin(alfa), radius1*cos(alfa))
    S2 = (dist-radius2*sin(beta), radius2*cos(beta))
    dx = S2[1] - S1[1]
    dy = S2[2] - S1[2]
    d = sqrt(dx^2 + dy^2)
    if (d^2 - (radius1 - radius2)^2) < 0
        return PATH_ERROR, nothing
    end
    l = sqrt(d^2 - (radius1 - radius2)^2)

    gamma = atan(dy/dx)
    xi = asin((radius1-radius2)/d)

    mu = mod2pi(gamma+xi + 3pi/2)
    if (S2[1] < S1[1])
        mu += pi
    end

    out = Vector{Float64}(undef,3)
    out[1] = mod2pi((mu - (alfa - pi/2))) * radius1
    out[2] = l
    out[3] = mod2pi((beta - pi/2 - mu)) * radius2

    return PATH_OK, out
end

function RSL_path(alfa, beta, dist, radius1, radius2)

    S1 = (radius1*sin(alfa), -radius1*cos(alfa))
    S2 = (dist-radius2*sin(beta), radius2*cos(beta))
    dx = S2[1] - S1[1]
    dy = S2[2] - S1[2]
    d = sqrt(dx^2 + dy^2)
    if (d^2 - (radius1 + radius2)^2) < 0
        return PATH_ERROR, nothing
    end
    l = sqrt(d^2 - (radius1 + radius2)^2)

    gamma = atan(dy/dx)
    xi = asin((radius1+radius2)/d)

    mu = mod2pi(gamma-xi + pi/2)

    out = Vector{Float64}(undef,3)
    out[1] = mod2pi(((alfa + pi/2) - mu)) * radius1
    out[2] = l
    out[3] = mod2pi(beta - pi/2 - (mu + pi)) * radius2

    return PATH_OK, out
end

function LSR_path(alfa, beta, dist, radius1, radius2)

    S1 = (-radius1*sin(alfa), radius1*cos(alfa))
    S2 = (dist+radius2*sin(beta), -radius2*cos(beta))

    dx = S2[1] - S1[1]
    dy = S2[2] - S1[2]
    d = sqrt(dx^2 + dy^2)
    if (d^2 - (radius1 + radius2)^2) < 0
        return PATH_ERROR, nothing
    end
    l = sqrt(d^2 - (radius1 + radius2)^2)

    gamma = atan(dy/dx)
    xi = asin((radius1+radius2)/d)

    mu = gamma + xi - pi/2

    out = Vector{Float64}(undef,3)
    out[1] = mod2pi(mu - (alfa - pi/2)) * radius1
    out[2] = l
    out[3] = mod2pi((mu + pi) - (beta + pi/2)) * radius2
    
    if isapprox(mod2pi((mu + pi) - (beta + pi/2)), 2pi, atol=1e-500)
        out[3] = 0.
    end

    return PATH_OK, out
end

"""
    get_configuration_in_turn(origin::Vector{Float64}, ang1::Float64, ang2::Float64, r::Float64, type::SegmentType)
Get vehicle configuration on a single arc.
Center of the arc is at coordinates [origin[1]-r*cos(ang1), origin[2]-r*sin(ang1)]
# Arguments
- `origin::Vector{Float64}`: vector of starting point: [x, y, _] or [x, y]
- `ang1::Float64`: starting angle in radians
- `ang2::Float64`: ending angle in radians
- `r::Float64`: radius of the arc
- `type::SegmentType`: wheater the arc is left/right turn
# Return
- `ret::Vector{Float64}`: vehicle configuration [x, y, theta]
"""
function get_configuration_in_turn(origin::Vector{Float64}, ang1::Float64, ang2::Float64, r::Float64, type::SegmentType)
    @assert type == L_SEG || type == R_SEG

    S = (origin[1]-r*cos(ang1), origin[2]-r*sin(ang1))
    x = S[1] + r * cos(ang2)
    y = S[2] + r * sin(ang2)
    theta = nothing
    if type == L_SEG
        theta = ang2 + pi/2
    else
        theta = ang2 - pi/2
    end
    ret::Vector{Float64} = [x, y, theta]
    return ret
end

"""
    get_configuration(path::DubinsPathR2, distance::Float64)
Get vehicle configuration at a given travelled distance
# Arguments
- `path::DubinsPathR2`: path on which the configuration is queried
- `distance::Float64`: distance from start at which the configuration is queried
# Return
- `ret::Vector{Float64}`: vehicle configuration [x, y, theta]
"""
function get_configuration(path::DubinsPathR2, distance::Float64)
    # Clip distance at end of the maneuver
    @assert distance >= 0. "Distance where configuration is queried must be greater than 0.0"
    if distance > sum(path.lengths)
        distance = sum(path.lengths)
    end

    if distance <= path.lengths[1]
        # Requested configuration is in the first segment
        if Integer(path.type) == Integer(LSR) || Integer(path.type) == Integer(LSL)
            ang1 = path.origin[3]-pi/2
            ang2 = ang1 + distance/path.r[1]
            r = path.r[1]
            x, y, theta = get_configuration_in_turn(path.origin, ang1, ang2, r, L_SEG)
        else
            ang1 = path.origin[3]+pi/2
            ang2 = ang1 - distance/path.r[1]
            r = path.r[1]
            x, y, theta = get_configuration_in_turn(path.origin, ang1, ang2, r, R_SEG)
        end
    elseif distance <= path.lengths[1] + path.lengths[2]
        # Requested conf is in second segment
        x = y = theta = nothing

        # Find end configuration of first segment
        if Integer(path.type) == Integer(LSR) || Integer(path.type) == Integer(LSL)
            ang1 = path.origin[3]-pi/2
            ang2 = ang1 + path.lengths[1]/path.r[1]
            r = path.r[1]
            x, y, theta = get_configuration_in_turn(path.origin, ang1, ang2, r, L_SEG)
        else
            ang1 = path.origin[3]+pi/2
            ang2 = ang1 - path.lengths[1]/path.r[1]
            r = path.r[1]
            x, y, theta = get_configuration_in_turn(path.origin, ang1, ang2, r, R_SEG)
        end

        # Find configuration within the second segment
        distance -= path.lengths[1]
        x, y = point_by_angle([x, y], theta, distance)
    else
        # Requested configuration is in the third segment
        x = y = theta= nothing

        # Find end configuration of the first segment
        if Integer(path.type) == Integer(LSR) || Integer(path.type) == Integer(LSL)
            ang1 = path.origin[3]-pi/2
            ang2 = ang1 + path.lengths[1]/path.r[1]
            r = path.r[1]
            x, y, theta = get_configuration_in_turn(path.origin, ang1, ang2, r, L_SEG)
        else
            ang1 = path.origin[3]+pi/2
            ang2 = ang1 - path.lengths[1]/path.r[1]
            r = path.r[1]
            x, y, theta = get_configuration_in_turn(path.origin, ang1, ang2, r, R_SEG)
        end

        # Find end configuration of the second segment
        x, y = point_by_angle([x, y], theta, path.lengths[2])
        distance -= (path.lengths[1] + path.lengths[2])

        # Find the configuration within the third segment
        if Integer(path.type) == Integer(LSL) || Integer(path.type) == Integer(RSL)
            ang1 = theta-pi/2
            ang2 = ang1 + distance/path.r[3]
            r = path.r[3]
            x, y, theta = get_configuration_in_turn([x, y], ang1, ang2, r, L_SEG)
        else
            ang1 = theta+pi/2
            ang2 = ang1 - distance/path.r[3]
            r = path.r[3]
            x, y, theta = get_configuration_in_turn([x, y], ang1, ang2, r, R_SEG)
        end
    end

    ret::Vector{Float64} = [x, y, theta]
    return ret
end 

###############################################################################
###################### SPEED PROFILE FUNCTIONS ################################
###############################################################################

"""
    speed_profile(final_path::Vector{DubinsPathR2}, params::Vector{Number}, speeds::Vector{Float64})
Sample single path for plotting

# Arguments
- `final_path::DubinsPathR2`: Vector of dubins paths
- `params::Vector{Float64}`: Vector of speed parameters: [velocity in minimal radius, maximal vehicle velocity, acceleration, deceleration]
- `speeds::Vector{Float64}` : Vector of speed values: [starting velocity, ending velocity]

# Return
- `times::Vector{Float64}`: times when the speed changes; Inf if invalid path
- `speed_at_time::Vector{Float64}`: speed at each time; Inf if invalid path
"""
function speed_profile(path::DubinsPathR2, params::Vector{Float64}, speeds::Vector{Float64})
    speed_profile(path, params[1], params[2], params[3], params[4], speeds[1], speeds[2])
end

function speed_profile(path::DubinsPathR2, v_min::Number, v_max::Number, a1::Float64, a2::Float64, v_i::Number, v_f::Number)
    # a1 = a_max, a2 = -a_min
    speeds = [v_i, min(v_max, speed_by_radius(path.r[2])), v_f]
    #Check if starting/ending speed is too high to make a turn of that radius
    if speed_by_radius(path.r[1]) < v_i || speed_by_radius(path.r[3]) < v_f
        return Inf, Inf
    end
    lengths = path.lengths

    times::Vector{Float64} = []
    speed_at_time::Vector{Float64} = []

    function push_vals(time, velocity)
        push!(times, time)
        push!(speed_at_time, velocity)
    end

    time = 0
    push_vals(time, speeds[1])

    time += lengths[1] / speeds[1]
    push_vals(time, speeds[1])
    # speedup in between
    sp_up = (speeds[2] - speeds[1])/a1
    sp_down = (speeds[2] - speeds[3])/a2
    sp_up_len = sp_up*(speeds[2] + speeds[1])/2
    sp_down_len = sp_down*(speeds[2] + speeds[3])/2

    if speeds[1] == speeds[3]
    #### CASE WITH SAME RADII

        if (sp_up_len + sp_down_len <= lengths[2])
        # full speedup (to maximum vs)
            time += sp_up
            push_vals(time, speeds[2])
            time += (lengths[2] - sp_down_len - sp_up_len) / speeds[2]
            push_vals(time, speeds[2])
            time += sp_down
            push_vals(time, speeds[3])
            time += lengths[3]/speeds[3]
            push_vals(time, speeds[3])

        else
        # small speedup (only to have enough time to speedup and speeddown)
        # find the intersect where it is necessary to change from speedup to speeddown
            vs = speeds[1] # start speed = end_speed
            s = lengths[2]
            # quadratic equation coeficients a, b, c, from the uniform speedup formula
            xa = a1 + a2 # xb = 0
            xc = -(vs^2)*xa - 2*a1*a2*s 
            det = -4*xa*xc
            if det < 0
                return Inf, Inf
            end

            vx = min(v_max, sqrt(det) / (2*xa)) # second root is smaller than beggining speed
            time += (vx - vs) / a1
            push_vals(time, vx)
            time += (vx - vs) / a2
            push_vals(time, vs)
            time += lengths[3] / vs
            push_vals(time, vs)
        end

    elseif speeds[1] < speeds[3]
    #### CASE WITH SPEEDUP TO THE LATTER SPEED
        # compute if straight part + latter part is long enough to speedup
        sp_up_finish = (speeds[3] - speeds[1])/a1
        sp_up_finish_len = sp_up_finish*(speeds[3] + speeds[1])/2

        # full speedup (to maximum vs)
        if (sp_up_len + sp_down_len) <= lengths[2]
            time += sp_up
            push_vals(time, speeds[2])
            time += (lengths[2] - sp_down_len - sp_up_len) / speeds[2]
            push_vals(time, speeds[2])
            time += sp_down
            push_vals(time, speeds[3])
            time += lengths[3] / speeds[3]
            push_vals(time, speeds[3])

        elseif sp_up_finish_len > lengths[2] && (sp_up_finish_len < lengths[2] + lengths[3])
            time += sp_up_finish
            push_vals(time, speeds[3])
            time += (lengths[3] + lengths[2] - sp_up_finish_len) / speeds[3]
            push_vals(time, speeds[3])

        elseif sp_up_finish_len < lengths[2] # small speedup (only to have enough time to speedup and speeddown)

            vs = speeds[1]
            ve = speeds[3]
            s = lengths[2]
            xa = a1 + a2 # xb = 0
            xc = -(vs^2)*a2 -(ve^2)*a1 - 2*a1*a2*s 
            det = -1 * (4*xa*xc) #
            if det < 0 
                @warn "det < 0"
                return Inf, Inf
            end

            vx = min(v_max, sqrt(det) / (2*xa)) # second root is smaller than previous speed - unnecessary
            if vx < ve
                return Inf, Inf
            end
            time += abs(vx - vs) / a1
            push_vals(time, vx)
            time += abs(vx - ve) / a2
            push_vals(time, ve)
            time += lengths[3] / ve
            push_vals(time, ve)
        else
            return Inf, Inf
        end

    elseif speeds[1] > speeds[3]
    #### CASE WITH SPEEDDOWN TO THE LATTER SPEED
        sp_down_finish = (speeds[1] - speeds[3])/a2
        sp_down_finish_len = sp_down_finish*(speeds[3] + speeds[1])/2

        if (sp_up_len + sp_down_len) < lengths[2]
            time += lengths[1] / speeds[1]
            push_vals(time, speeds[1])
            time += sp_up
            push_vals(time, speeds[2])
            time += (lengths[2] - sp_down_len - sp_up_len) / speeds[2]
            push_vals(time, speeds[2])
            time += sp_down
            push_vals(time, speeds[3])
            time += lengths[3] / speeds[3]
            push_vals(time, speeds[3])
        elseif sp_down_finish_len > lengths[2] && sp_down_finish_len < (lengths[1] + lengths[2])
            time -= lengths[1] / speeds[1]
            pop!(times)
            pop!(speed_at_time)
            time += (lengths[1] + lengths[2] - sp_down_finish_len) / speeds[1]
            push_vals(time, speeds[1])
            time += sp_down_finish
            push_vals(time, speeds[3])
            time += lengths[3] / speeds[3]
            push_vals(time, speeds[3])
        elseif sp_down_finish_len < lengths[2]
            vs = speeds[1]
            ve = speeds[3]
            s = lengths[2]
            xa = a1 + a2 # xb = 0
            xc = -(vs^2)*a2 -(ve^2)*a1 - 2*a1*a2*s 
            det = -1 * (4*xa*xc)
            if det < 0
                return Inf, Inf
            end
            vx = min(v_max, sqrt(det) / (2*xa)) # second root is smaller than previous speed - unnecessary
            if vx < ve
                return Inf, Inf
            end

            time += lengths[1] / speeds[1]
            push_vals(time, speeds[1])
            time += abs(vx - vs) / a1
            push_vals(time, vx)
            time += abs(vx - ve) / a2
            push_vals(time, ve)
            time += lengths[3] / speeds[3]
            push_vals(time, speeds[3])
        else
            return Inf, Inf
        end

    end
    return times, speed_at_time
end

end