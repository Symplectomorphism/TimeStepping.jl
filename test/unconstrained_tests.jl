using Test

include("bouncing_ball.jl")

@testset "Penetration Tests" begin
    @test all(getindex.(bouncing_ball.q,2) .>= -1e-3)
end

@testset "Energy Tests" begin
    energy(q, u) = 1/2 * mass * (u[1]^2 + u[2]^2) + mass * g * q[2]
    initial_energy = energy(bouncing_ball.q[1], bouncing_ball.u[1])

    @test all(energy.(bouncing_ball.q, bouncing_ball.u) .<= initial_energy)
end