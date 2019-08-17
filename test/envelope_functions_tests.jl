using DataFrames
using NamedArrays
using SchumakerSpline
using Random

bounds = (0.0, 1.0)
x_name = :alpha
evaluations = DataFrame()

function func_inc1(x)
    feval = NamedArray([x, -sqrt(x)])
    setnames!(feval, ["utility", "x_E"], 1)
    return feval
end

function func_inc2(x)
    feval = NamedArray([sqrt(x), sqrt(x)])
    setnames!(feval, ["utility", "x_E"], 1)
    return feval
end

function func_inc3(x)
    feval = NamedArray([3*(x^2), sqrt(x)])
    setnames!(feval, ["utility", "x_E"], 1)
    return feval
end

function dec1(x)
    feval = NamedArray([-sqrt(x), sqrt(x)])
    setnames!(feval, ["utility", "x_E"], 1)
    return feval
end
function dec2(x)
    feval = NamedArray([-3*x^2 + 1, sqrt(x)])
    setnames!(feval, ["utility", "x_E"], 1)
    return feval
end

const_Array = NamedArray([0.1, sqrt(0.2)])
setnames!(const_Array, ["utility", "x_E"], 1)

# Testing the construction of a dataframe of functions
dd_function = make_function_dataframe([:func_inc1, :func_inc2, :func_inc3, :dec1, :dec2, :constant],
                                      [func_inc1, func_inc2, func_inc3, dec1, dec2, const_Array],
                                      [:Increasing, :Increasing, :Increasing, :Decreasing, :Decreasing, :Constant])
grid                = reverse(0.0001:0.01:0.9999)
output_for_envelope = :utility
x_name              = :alpha
env                 = Envelope(dd_function, grid, output_for_envelope,x_name)

x = 0.34

evaluate(env, x)
evaluate(env, x; derivative = 1)
evaluate(env, x, derivative = 2)
evaluate(env, x, derivative = 3)

evaluate_integral(env, 0.2, 0.9)

sample_envelope_given_func(env, :func_inc3, :utility, 0.9)

twister = MersenneTwister(1)
var = :utility
func_name, x, y = sample_envelope(env, var, twister; derivative = 0)
x < 1
