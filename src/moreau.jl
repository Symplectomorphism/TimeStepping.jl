mutable struct Moreau
    dynamics::Function
    gap::Function                   # Gap function -- provide the force matrices
    M::Array{Float64, 2}
    h::Array{Float64, 1}
    tA::Float64
    tM::Float64
    tE::Float64
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

function Moreau(n::Int, gap::Function, dynamics::Function, q::Vector, u::Vector)
    Δt = 1e-3
    qA = q
    uA = u
    qM = _compute_mid_displacements(qA, uA, Δt)
    M, h = dynamics(qM, uA)
    qE = zeros(Float64, n)
    uE = zeros(Float64, n)
    tA, tM, tE = zeros(Float64, 3)
    tM = _compute_mid_time(tA, Δt)
    tE = _compute_mid_time(tM, Δt)
    ε = 0.5
    g, W = gap(qA, uA)
    H = zeros(Int, 1)
    Λ = zeros(Float64, length(H))

    Moreau(dynamics, gap, M, h, tA, tM, tE, qA, uA, qM, qE, uE, g, W, Λ, H, Δt, ε)
end


function _compute_mid_time(t::Number, Δt::Number)
    return t + 1/2*Δt
end

function _compute_mid_displacements(q::Vector, u::Vector, Δt::Number)
    return q + 1/2*Δt*u
end