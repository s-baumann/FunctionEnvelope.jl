using FunctionEnvelope
using Test

# Run tests
println("Basic Tests")
@time @test include("basic_tests.jl")
println("Full Tests")
@time @test include("envelope_functions_tests.jl")
println("Basic Tests")
@time @test include("one_function_test.jl")
