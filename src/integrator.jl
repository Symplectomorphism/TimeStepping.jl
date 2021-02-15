mutable struct Integrator
    t::Array{Float64, 1}
    q::Array{Array{Float64, 1}, 1}
    u::Array{Array{Float64, 1}, 1}
    Λ::Array{Array{Float64, 1}, 1}
    dynamics::Function
    gap::Function
    m::Moreau
    Δt::Float64
end

function Integrator(gap::Function, dynamics::Function, q0::Vector, u0::Vector; Δt::Float64=1e-3)
    t = Array{Float64, 1}()
    q = Array{Array{Float64, 1},1}()
    u = Array{Array{Float64, 1},1}()
    Λ = Array{Array{Float64, 1},1}()
    push!(t, 0.0)
    push!(q, q0)
    push!(u, u0)
    m = Moreau(gap, dynamics, q[1], u[1], Δt)

    Integrator(t, q, u, Λ, dynamics, gap, m, Δt)
end


function integrate(system::Integrator, final_time::Float64)
    for time in range(system.t[1]; step=system.Δt, stop=final_time)
        push!(system.t, time)
        step(system.m)
        push!(system.q, system.m.qE)
        push!(system.u, system.m.uE)
        push!(system.Λ, system.m.Λ)
        set_state(system.m, system.q[end], system.u[end])
    end
end