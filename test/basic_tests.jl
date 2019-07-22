using FunctionEnvelope
using DataFrames
using NamedArrays
using SchumakerSpline

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

# Testing basic evaluations
dd = evaluate_function(func_inc1,0.3,:func_inc1, :example_x)
size(dd) == (1, 4)
output, dd2 = evaluate_function_and_add_to_dd(dd, func_inc2,0.4,:func_inc2, :example_x, :utility)
size(dd2) == (2, 4)
output, dd3 = evaluate_function_and_add_to_dd(dd2, dec1,0.4,:dec1, :example_x, :utility)
size(dd3) == (3, 4)
output, dd4 = evaluate_function_and_add_to_dd(dd3, const_Array,0.4,:const_Array, :example_x, :utility)
size(dd4) == (4, 4)

# Testing the construction of a dataframe of functions
dd_function = make_function_dataframe([:func_inc1, :func_inc2, :func_inc3, :dec1, :dec2, :constant],
                                      [func_inc1, func_inc2, func_inc3, dec1, dec2, const_Array],
                                      [:Increasing, :Increasing, :Increasing, :Decreasing, :Decreasing, :Constant])
dd_function[2,:func] == func_inc2

# Testing the crawler.
step = 0.11
limits = (0.0001, 0.9999)
step_right_to_left = true
x_name = :alpha
output_for_envelope = :utility

# Getting Initial evaluations
dd_eval = initial_evaluations(dd_function, limits, step_right_to_left, x_name, output_for_envelope)

# Testing
target = 0.8
max_vals =  maximum_possible_value(target, dd_function, dd_eval, output_for_envelope, x_name)

# Initial evaluations
grid = reverse(0.0001:step:0.9999)
dd_eval = crawler(dd_function, grid, x_name, output_for_envelope; do_all_at_ends = true)


x_array = sort(unique(dd_eval[x_name]))
bracketing_parameter = 0.8
max_interval_width = 0.01
func_name = :func_inc2
interval = Tuple{Float64,Float64}([0.01, 0.8])
type = :Increasing
constant = 0.2
dd_eval_rows = size(dd_eval)[1]
# This should have one intercept at x = 0.2^2 = 0.04
interc, dd_eval = get_intercept_with_constant(dd_eval, func_name, type, constant, dd_function, interval, x_name, output_for_envelope, bracketing_parameter, max_interval_width)
size(dd_eval)[1] > dd_eval_rows
abs(interc - 0.04) < max_interval_width
# This should have one intercept at x = 0.2
func_name = :func_inc1
interc, dd_eval = get_intercept_with_constant(dd_eval, func_name, type, constant, dd_function, interval, x_name, output_for_envelope, bracketing_parameter, max_interval_width)
abs(interc - 0.2) < max_interval_width

# These should intercept at x = 1/3
increasing_function = :func_inc3
decreasing_function = :func_inc1
interval = Tuple{Float64,Float64}([0.1,0.6])
interp, dd_eval = get_intercept(dd_eval, increasing_function, decreasing_function, dd_function, interval,
                                x_name, output_for_envelope, bracketing_parameter, max_interval_width)
abs(interp - (1/3)) < max_interval_width

# These should intercept at x = (-1 + sqrt(13))/(6)
increasing_function = :func_inc1
decreasing_function = :dec2
interval = Tuple{Float64,Float64}([0.1,0.6])
interp, dd_eval = get_intercept(dd_eval, increasing_function, decreasing_function, dd_function, interval,
                                x_name, output_for_envelope, bracketing_parameter, max_interval_width)
analytical_intercept = (-1 + sqrt(13))/(6) # Just do quadratic formula to get this.
abs(interp - analytical_intercept) < max_interval_width



# Based on the above functions
dd_interval, dd_eval = get_upper_envelope(dd_eval, dd_function, x_name, output_for_envelope, bracketing_parameter, max_interval_width; grid = reverse(collect(grid)))
env = Envelope(dd_interval, dd_eval, dd_function, x_name, output_for_envelope)


plt = plot(env; points = true, legend = :none)
plt
plot(plt; legend=:best)
