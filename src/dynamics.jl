mutable struct Dynamics
    M::Array{Float64, 2}
    h::Array{Float64, 1}
    Δt::Float64
    qA::Array{Float64, 1}
    uA::Array{Float64, 1}
    qM::Array{Float64, 1}
    qE::Array{Float64, 1}
    uE::Array{Float64, 1}
    W::Array{Float64, 2}            # Jacobian transpose in the normal direction
    Λ::Array{Float64, 1}            # Normal contact force
    g::Function                     # Gap function -- provide the force matrices
    ε::Float64                      # Coefficient of restitution (normal dir.)
end


function Dynamics(n::Int)
    M = Matrix(1.0I, n, n)
    h = zeros(Float64, n)
    Δt = 1e-3
end