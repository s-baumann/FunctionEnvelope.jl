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
const_Array = NamedArray([0.1, sqrt(0.2)])
setnames!(const_Array, ["utility", "x_E"], 1)

# Testing the construction of a dataframe of functions
dd_function = make_function_dataframe([:func_inc1], Array{Any,1}([func_inc1]), [:Increasing])
dd_function2 = make_function_dataframe([:const_Array], Array{Any,1}([const_Array]), [:Constant])

#
grid                = reverse(0.0001:0.01:0.9999)
output_for_envelope = :utility
x_name              = :alpha
env                 = Envelope(dd_function, grid, output_for_envelope,x_name)
isa(env, Envelope)
env2                = Envelope(dd_function2, grid, output_for_envelope,x_name)
isa(env2, Envelope)
