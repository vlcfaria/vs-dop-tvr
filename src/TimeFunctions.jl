module TimeFunctions

function sine_wave(a::Number, b::Number, f::Number, d::Number)
    return t -> (a + b * sin(2 * pi * f * t + d))
end

end