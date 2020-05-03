module FunctionEnvelope

using BSON
using DataFrames
using NamedArrays
using SchumakerSpline
using Plots
using Distributions
using Random


include("1_structs.jl")
export make_function_dataframe
include("2_evaluate_functions.jl")
export evaluate_function, evaluate_function_and_add_to_dd, get_function_value
include("3_crawler.jl")
export initial_evaluations, maximum_possible_value, get_func, top_function_at_point, crawler
export get_intercept_with_constant, get_intercept, get_upper_envelope
include("4_envelope.jl")
export Envelope, evaluate, evaluate_integral, sample_envelope_given_func, sample_envelope
include("5_plotting.jl")
export plot
end
