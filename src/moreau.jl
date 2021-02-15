mutable struct Moreau
    dynamics::Function
    gap::Function                   # Gap function -- provide the force matrices
    M::Array{Float64, 2}
    h::Array{Float64, 1}
    qA::Array{Float64, 1}
    uA::Array{Float64, 1}
    qM::Array{Float64, 1}
    qE::Array{Float64, 1}
    uE::Array{Float64, 1}
    g::Array{Float64, 1}            # Normal distance of contact
    W::Array{Float64, 2}            # Jacobian transpose in the normal direction
    Λ::Array{Float64, 1}            # Normal contact force
    H::Array{Int, 1}                # Which contacts are active?
    Δt::Float64
    ε::Float64                      # Coefficient of restitution (normal dir.)
end

"""
n: degrees of freedom of the system
m: number of unilateral constraints
"""

function Moreau(n::Int, gap::Function, dynamics::Function)
    qA = zeros(Float64, n)
    uA = zeros(Float64, n)
    qM = zeros(Float64, n)
    qE = zeros(Float64, n)
    uE = zeros(Float64, n)
    M, h = dynamics(qA, uA)
    Δt = 1e-3
    ε = 0.5
    g, W = gap(qA, uA)
    H = zeros(Int, 0)
    Λ = zeros(Float64, length(H))

    Moreau(dynamics, gap, M, h, qA, uA, qM, qE, uE, g, W, Λ, H, Δt, ε)
end

