using SafeTestsets

@safetestset "Unconstrained tests" begin include("unconstrained_tests.jl") end
@safetestset "Electromechanical tests" begin include("electromechanical_tests.jl") end
@safetestset "Constrained tests" begin include("constrained_tests.jl") end