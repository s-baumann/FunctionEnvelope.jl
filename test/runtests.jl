using FunctionEnvelope
using Test

# Run tests
println("Basic Tests")
@time @test include("basic_tests.jl")
