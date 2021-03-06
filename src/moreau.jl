const MOI = MathOptInterface
const SUCCESS = [MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED]

mutable struct Moreau{SpecializedType}
    dynamics::Function
    gap::Function                   # Gap function -- provide the force matrices
    hcon::Function                  # Holonomic constraint function
    jac::Function                   # Jacobian of hcon: del hcon / del q
    jacdot::Function                # Time-derivative of jac
    M::Array{SpecializedType, 2}
    h::Array{SpecializedType, 1}
    ϕ::Array{SpecializedType, 1}
    J::Array{SpecializedType, 2}
    Jdot::Array{SpecializedType, 2}
    tA::SpecializedType
    tM::SpecializedType
    tE::SpecializedType
    qA::Array{SpecializedType, 1}
    uA::Array{SpecializedType, 1}
    qM::Array{SpecializedType, 1}
    qE::Array{SpecializedType, 1}
    uE::Array{SpecializedType, 1}
    extra_state::Array{SpecializedType, 1}  
    g::Array{SpecializedType, 1}            # Normal distance of contact
    W::Array{SpecializedType, 2}            # Jacobian transpose in the normal direction (changes during integration)
    Wn::Array{SpecializedType, 2}           # Jacobian transpose in the normal direction (fixed, does not change)
    Λ::Array{SpecializedType, 1}            # Normal contact force
    μ::Array{SpecializedType, 1}            # Holonomic constraint forces
    H::SortedSet{Int64, Base.Order.ForwardOrdering} # Which contacts are active?
    Δt::SpecializedType
    ε::SpecializedType                      # Coefficient of restitution (normal dir.)
    ϵ::SpecializedType                      # Baumgarte stabilization constant
end

"""
n: degrees of freedom of the system
m: number of unilateral constraints
"""

function Moreau(gap::Function, dynamics::Function, q::AbstractArray, u::AbstractArray, 
        Δt::Wildcard=1e-3, ε::Wildcard=0.5) where {Wildcard <: Real}
    qA = q
    uA = u
    extra_state = Wildcard[]
    n = length(qA)
    qM = _compute_mid_displacements(qA, uA, Δt)
    M, h = dynamics(qM, uA)
    ϕ = Wildcard[]
    J = Array{Wildcard, 2}(undef, 0, 0)
    Jdot = Array{Wildcard, 2}(undef, 0, 0)
    qE = zeros(Wildcard, n)
    uE = zeros(Wildcard, n)
    tA, tM, tE = zeros(Wildcard, 3)
    tM = _compute_mid_time(tA, Δt)
    tE = _compute_mid_time(tM, Δt)
    ϵ = 0.01
    g, Wn = gap(qA, uA)
    W = Wn
    H = SortedSet(Int[])
    Λ = zeros(Wildcard, length(H))
    μ = zeros(Wildcard, length(ϕ))

    Moreau{Wildcard}(dynamics, gap, x->x, x->x, x->x,
        M, h, ϕ, J, Jdot, tA, tM, tE, qA, uA, qM, qE, uE, extra_state, 
        g, W, Wn, Λ, μ, H, Δt, ε, ϵ)
end

function Moreau(gap::Function, dynamics::Function, hcon::Function, 
        jac::Function, jacdot::Function, q::AbstractArray, u::AbstractArray, 
        Δt::Wildcard=1e-3, ε::Wildcard=0.5) where {Wildcard <: Real}
    qA = q
    uA = u
    extra_state = Wildcard[]
    n = length(qA)
    qM = _compute_mid_displacements(qA, uA, Δt)
    M, h = dynamics(qM, uA)
    ϕ = hcon(qM)
    J = jac(qM)
    Jdot = jacdot(qM, uA)
    qE = zeros(Wildcard, n)
    uE = zeros(Wildcard, n)
    tA, tM, tE = zeros(Wildcard, 3)
    tM = _compute_mid_time(tA, Δt)
    tE = _compute_mid_time(tM, Δt)
    ϵ = 0.01
    g, Wn = gap(qA, uA)
    W = Wn
    H = SortedSet(Int[])
    Λ = zeros(Wildcard, length(H))
    μ = zeros(Wildcard, length(ϕ))

    Moreau{Wildcard}(dynamics, gap, hcon, jac, jacdot,
        M, h, ϕ, J, Jdot, tA, tM, tE, qA, uA, qM, qE, uE, extra_state, 
        g, W, Wn, Λ, μ, H, Δt, ε, ϵ)
end


function Moreau(gap::Function, dynamics::Function, q::AbstractArray, u::AbstractArray, 
    extra_state::AbstractArray, Δt::Wildcard=1e-3, ε::Wildcard=0.1) where {Wildcard <: Real}
    qA = q
    uA = u
    n = length(qA)
    qM = _compute_mid_displacements(qA, uA, Δt)
    M, h = dynamics(qM, uA, extra_state)
    ϕ = Wildcard[]
    J = Array{Wildcard, 2}(undef, 0, 0)
    Jdot = Array{Wildcard, 2}(undef, 0, 0)
    qE = zeros(Wildcard, n)
    uE = zeros(Wildcard, n)
    tA, tM, tE = zeros(Wildcard, 3)
    tM = _compute_mid_time(tA, Δt)
    tE = _compute_mid_time(tM, Δt)
    ϵ = 0.01
    g, Wn = gap(qA, uA)
    W = Wn
    H = SortedSet(Int[])
    Λ = zeros(Wildcard, length(H))
    μ = zeros(Wildcard, length(ϕ))

    Moreau{Wildcard}(dynamics, gap, x->x, x->x, x->x,
        M, h, ϕ, J, Jdot, tA, tM, tE, qA, uA, qM, qE, uE, extra_state, 
        g, W, Wn, Λ, μ, H, Δt, ε, ϵ)
end


function _compute_mid_time(t::Number, Δt::Real)
    return t + 1/2*Δt
end

function _compute_mid_displacements(q::AbstractArray, u::AbstractArray, Δt::Real)
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
    m.W = zeros(eltype(m.qA), n, 0)
    for i in m.H
        m.W = hcat(m.W, m.Wn[:,i])
    end
    m.Λ = zeros(eltype(m.qA), length(m.g))
end

function step(m::Moreau)
    if isempty(m.ϕ)
        step_unconstrained(m)
    else
        step_constrained(m)
    end
end

function step_unconstrained(m::Moreau)
    _compute_index_set(m)
    Minv = inv(m.M)

    model = Model(Mosek.Optimizer)
    set_optimizer_attribute(model, "QUIET", true)
    set_optimizer_attribute(model, "INTPNT_CO_TOL_DFEAS", 1e-7)
    @variable(model, q[1:length(m.qA)])
    @variable(model, u[1:length(m.uA)])
    @variable(model, λ[1:length(m.H)] >= 0)

    if length(m.W) != 0
        """ The following does not work because of positive semidefiniteness...
        ξ = m.W' * ( u + m.ε * m.uA )
        @constraint(model, ξ .>= 0)
        @objective(model, Min, dot(λ, ξ))
        """
        A = m.W' * Minv * m.W
        b = m.W' * Minv * m.h * m.Δt + (1+m.ε) * m.W' * m.uA
        @constraint(model, b .+ A*λ .>= 0)
        @constraint(model, m.M * (u - m.uA) .== m.W * λ + m.h * m.Δt)
        @objective(model, Min, dot(λ, b .+ A*λ))
    else
        @constraint(model, λ .== 0)
        @constraint(model, m.M * (u - m.uA) .== m.h * m.Δt)
        @objective(model, Min, 0)
    end
    @constraint(model, q .== _compute_mid_displacements(m.qM, u, m.Δt))
    optimize!(model)
    if !any( JuMP.termination_status(model) .== SUCCESS )
        @warn "termination status: $(JuMP.termination_status(model))."
    else
        m.qE = JuMP.value.(q)
        m.uE = JuMP.value.(u)
        if length(m.W) != 0
            i = 1
            for k in 1:length(m.g)
                if k in m.H
                    m.Λ[k] = JuMP.value.(λ)[i]
                    i += 1
                else
                    m.Λ[k] = 0.0
                end
            end
        else
            m.Λ = zeros(eltype(m.qA), length(m.g))
        end
    end
end

function step_constrained(m::Moreau)
    _compute_index_set(m)
    Minv = inv(m.M)
    hJ = -m.Jdot * m.uA - 2/m.ϵ*m.J*m.uA - 1/m.ϵ/m.ϵ*m.ϕ
    Mhat = inv( m.J * Minv * m.J' )

    model = Model(Mosek.Optimizer)
    set_optimizer_attribute(model, "QUIET", true)
    set_optimizer_attribute(model, "INTPNT_CO_TOL_DFEAS", 1e-7)
    @variable(model, q[1:length(m.qA)])
    @variable(model, u[1:length(m.uA)])
    @variable(model, λ[1:length(m.H)] >= 0)

    if length(m.W) != 0
        A = m.W' * Minv * ( I - m.J' * Mhat * m.J * Minv ) * m.W
        b = m.W' * Minv * ( m.h + m.J' * Mhat *(hJ - m.J * Minv * m.h) ) * m.Δt + (1+m.ε) * m.W' * m.uA
        μ = Mhat * ( hJ - m.J * Minv * m.h - 1/m.Δt * m.J * Minv * m.W * λ )
        @constraint(model, b .+ A*λ .>= 0)
        @constraint(model, m.M * (u - m.uA) .== m.W * λ + (m.h + m.J' * μ) * m.Δt)
        @objective(model, Min, dot(λ, b .+ A*λ))
    else
        μ = Mhat * ( hJ - m.J * Minv * m.h )
        @constraint(model, λ .== 0)
        @constraint(model, m.M * (u - m.uA) .== (m.h + m.J' * μ) * m.Δt)
        @objective(model, Min, 0)
    end
    @constraint(model, m.J*(u - m.uA) .== hJ * m.Δt)
    @constraint(model, q .== _compute_mid_displacements(m.qM, u, m.Δt))
    optimize!(model)
    if !any( JuMP.termination_status(model) .== SUCCESS )
        @warn "termination status: $(JuMP.termination_status(model))."
    else
        m.qE = JuMP.value.(q)
        m.uE = JuMP.value.(u)
        if length(m.W) != 0
            i = 1
            for k in 1:length(m.g)
                if k in m.H
                    m.Λ[k] = JuMP.value.(λ)[i]
                    i += 1
                else
                    m.Λ[k] = 0.0
                end
            end
            m.μ = Mhat * ( hJ - m.J * Minv * m.h - 1/m.Δt * m.J * Minv * m.W * m.Λ )
        else
            m.Λ = zeros(eltype(m.qA), length(m.g))
            m.μ = Mhat * ( hJ - m.J * Minv * m.h )
        end
    end
end

function set_state(m::Moreau, q::AbstractArray, u::AbstractArray)
    m.qA = q
    m.uA = u
    n = length(m.qA)
    m.qM = _compute_mid_displacements(m.qA, m.uA, m.Δt)
    m.M, m.h = m.dynamics(m.qM, m.uA)
    m.g, m.W = m.gap(m.qM, m.uA)
    if !isempty(m.ϕ)
        m.ϕ = m.hcon(m.qM)
        m.J = m.jac(m.qM)
        m.Jdot = m.jacdot(m.qM, m.uA)
    end
    _compute_index_set(m)
end

function set_state(m::Moreau, q::AbstractArray, u::AbstractArray, extra_state::AbstractArray)
    m.qA = q
    m.uA = u
    m.extra_state = extra_state
    n = length(m.qA)
    m.qM = _compute_mid_displacements(m.qA, m.uA, m.Δt)
    m.M, m.h = m.dynamics(m.qM, m.uA, m.extra_state)
    m.g, m.W = m.gap(m.qM, m.uA)
    if !isempty(m.ϕ)
        m.ϕ = m.hcon(m.qM)
        m.J = m.jac(m.qM)
        m.Jdot = m.jacdot(m.qM, m.uA)
    end
    _compute_index_set(m)
end