struct Integrator{T}
    t::Array{T, 1}
    q::Array{Array{T, 1}, 1}
    u::Array{Array{T, 1}, 1}
    Λ::Array{Array{T, 1}, 1}
    m::Moreau
    Δt::T
end


function Integrator(gap::Function, dynamics::Function, q0::AbstractArray, 
        u0::AbstractArray; Δt::S=1e-3*1f0) where {S <: Real}
    t = Array{S, 1}()
    q = Array{Array{S, 1},1}()
    u = Array{Array{S, 1},1}()
    Λ = Array{Array{S, 1},1}()
    push!(t, 0.0)
    push!(q, q0)
    push!(u, u0)
    m = Moreau(gap, dynamics, q[1], u[1], Δt)
    set_state(m, q0, u0)

    Integrator{S}(t, q, u, Λ, m, Δt)
end


function Integrator(gap::Function, dynamics::Function, hcon::Function, 
        jac::Function, jacdot::Function, q0::AbstractArray, u0::AbstractArray; Δt::S=1e-3) where {S <: Real}
    t = Array{S, 1}()
    q = Array{Array{S, 1},1}()
    u = Array{Array{S, 1},1}()
    Λ = Array{Array{S, 1},1}()
    push!(t, 0.0)
    push!(q, q0)
    push!(u, u0)
    m = Moreau(gap, dynamics, hcon, jac, jacdot, q[1], u[1], Δt)
    set_state(m, q0, u0)

    Integrator{S}(t, q, u, Λ, m, Δt)
end


function integrate(system::Integrator, final_time::S) where {S <: Real}
    for time in range(system.t[1]; step=system.Δt, stop=final_time)
        push!(system.t, time)
        step(system.m)
        push!(system.q, system.m.qE)
        push!(system.u, system.m.uE)
        push!(system.Λ, system.m.Λ)
        set_state(system.m, system.q[end], system.u[end])
    end
end