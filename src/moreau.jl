const MOI = MathOptInterface
const SUCCESS = [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED]

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
    H::SortedSet{Int64, Base.Order.ForwardOrdering} # Which contacts are active?
    Δt::Float64
    ε::Float64                      # Coefficient of restitution (normal dir.)
end

"""
n: degrees of freedom of the system
m: number of unilateral constraints
"""

function Moreau(gap::Function, dynamics::Function, q::Vector, u::Vector, Δt::Float64=1e-3)
    qA = q
    uA = u
    n = length(qA)
    qM = _compute_mid_displacements(qA, uA, Δt)
    M, h = dynamics(qM, uA)
    qE = zeros(Float64, n)
    uE = zeros(Float64, n)
    tA, tM, tE = zeros(Float64, 3)
    tM = _compute_mid_time(tA, Δt)
    tE = _compute_mid_time(tM, Δt)
    ε = 0.5
    g, W = gap(qA, uA)
    H = SortedSet(Int[])
    Λ = zeros(Float64, length(H))

    Moreau(dynamics, gap, M, h, tA, tM, tE, qA, uA, qM, qE, uE, g, W, Λ, H, Δt, ε)
end


function _compute_mid_time(t::Number, Δt::Number)
    return t + 1/2*Δt
end

function _compute_mid_displacements(q::Vector, u::Vector, Δt::Number)
    return q + 1/2*Δt*u
end

function _compute_index_set(m::Moreau)
    contact_threshold = 1e-4
    i=1
    for g in m.g
        if g <= contact_threshold
            push!(m.H, i)
        else
            try
                pop!(m.H, i)
            catch e
            end
        end
        i += 1
    end
    _update_force_matrix(m)
end

function _update_force_matrix(m::Moreau)
    n = length(m.qA)
    temp = zeros(Float64, n, 0)
    for i in m.H
        temp = hcat(temp, m.W[:,1])
    end
    m.W = temp
    m.Λ = zeros(Float64, length(m.g))
end

function _solve_LCP(m::Moreau)
    Minv = inv(m.M)
    if length(m.W) != 0
        A = m.W' * Minv * m.W
        b = m.W' * Minv * m.h * m.Δt + (1+m.ε) * m.W' * m.uA

        model = Model(PATHSolver.Optimizer)
        set_optimizer_attribute(model, "output", "no")
        @variable(model, x[1:length(m.Λ)] >= 0)
        @constraint(model, A*x .+ b ⟂ x)
        optimize!(model)
        if !any( JuMP.termination_status(model) .== SUCCESS )
            @warn "termination status: $(JuMP.termination_status(model))."
        else
            m.Λ = JuMP.value.(x)
        end
    else
        m.Λ = zeros(Float64, length(m.g))
    end
    return Minv
end

function step(m::Moreau)
    _compute_index_set(m)
    Minv = _solve_LCP(m)
    if length(m.W) != 0
        m.uE = Minv * m.W * m.Λ + Minv * m.h * m.Δt + m.uA
    else
        m.uE = Minv * m.h * m.Δt + m.uA
    end
    m.qE = _compute_mid_displacements(m.qM, m.uE, m.Δt)
end

function step2(m::Moreau)
    _compute_index_set(m)
    Minv = inv(m.M)

    model = Model(Mosek.Optimizer)
    set_optimizer_attribute(model, "QUIET", true)
    set_optimizer_attribute(model, "INTPNT_CO_TOL_DFEAS", 1e-7)
    @variable(model, q[1:length(m.qA)])
    @variable(model, u[1:length(m.uA)])
    @variable(model, λ[1:length(m.Λ)] >= 0)

    if length(m.W) != 0
        A = m.W' * Minv * m.W
        b = m.W' * Minv * m.h * m.Δt + (1+m.ε) * m.W' * m.uA
        @constraint(model, A*λ + b .>= 0)
        @constraint(model, m.M * (u - m.uA) .== m.W * λ + m.h * m.Δt)
        @objective(model, Min, dot(λ, b .+ A*λ))
    else
        @constraint(model, λ .== 0)
        @constraint(model, m.M * (u - m.uA) .== m.h * m.Δt)
        @objective(model, Min, 0)
    end
    @constraint(model, q .== m.qM + 1/2 * m.Δt * u)
    optimize!(model)
    if !any( JuMP.termination_status(model) .== SUCCESS )
        @warn "termination status: $(JuMP.termination_status(model))."
    else
        m.qE = JuMP.value.(q)
        m.uE = JuMP.value.(u)
        if length(m.W) != 0
            m.Λ = JuMP.value.(λ)
        else
            m.Λ = zeros(Float64, length(m.g))
        end
    end
end

function set_state(m::Moreau, q::Vector, u::Vector)
    m.qA = q
    m.uA = u
    n = length(m.qA)
    m.qM = _compute_mid_displacements(m.qA, m.uA, m.Δt)
    m.M, m.h = m.dynamics(m.qM, m.uA)
    m.g, m.W = m.gap(m.qA, m.uA)
    _compute_index_set(m)
end