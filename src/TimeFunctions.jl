module TimeFunctions

using FunctionWrappers
import FunctionWrappers: FunctionWrapper

const TimeFunc = FunctionWrapper{Float64, Tuple{Float64}}

function sine_wave(a::Number, b::Number, f::Number, d::Number)
    if a - abs(b) < 0
        println("Warning: sine wave with negative values")
    end

    return t -> (a + b * sin(2 * pi * f * t + d))
end

end