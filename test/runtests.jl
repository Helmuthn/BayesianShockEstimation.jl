using BayesianShockEstimation
using Test
using SafeTestsets

@testset "BayesianShockEstimation.jl" begin
    @safetestset "rts.jl" begin
        include("rts_tests.jl")
    end
    @safetestset "burgers.jl" begin
        include("burgers_tests.jl")
    end
end
