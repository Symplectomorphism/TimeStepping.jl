"""
Magnetic levitation system
q = (z), u = (z1dot), extra_state = i
g(q) = (q-0.05, 0.20-q)
W = del g / del q = [1 -1]
M = m*Matrix(1.0f0*I, 1, 1)
h = m*g*ones(Float32, 1)
"""

using LinearAlgebra
using PyPlot
using TimeStepping

par = Dict{Symbol, Float32}(
    :m=>0.01187, :g=>9.81, :C=>1.4e-4, :L1=>0.65, :R=>28.7, :qd=>0.1, :tf=>4.0,
    :Δt=>1e-3, :ε=>0.1
)

δ = convert.(Float32, [0.05, 0.20])       # Hard stops of the maglev system.
q0 = convert.(Float32, [δ[1]])            # Initial ball position
u0 = convert.(Float32, [0.0])             # Initial ball velocity
extra_state0 = convert.(Float32, [0.0])   # Initial magnet current    (5.784 to lift-off)

function gap(q::AbstractArray, u::AbstractArray)
    return [q[1].-δ[1], δ[2].-q[1]], convert.(Float32, [1.0 -1])
end

function dynamics(q::AbstractArray, u::AbstractArray, extra_state::AbstractArray)
    return par[:m]*Matrix(1.0f0*I, 1, 1), 
        (par[:m]*par[:g] .- par[:C]*(extra_state[1]./q[1]).^2)*ones(Float32, 1)
end

function fl_controller(extra_state::AbstractArray, q::AbstractArray, u::AbstractArray)
    m = par[:m]
    g = par[:g]
    C = par[:C]
    L1 = par[:L1]
    R = par[:R]
    qd = par[:qd]

    x = (q[1], u[1], extra_state[1])
    ξ = (x[1]-qd, x[2], g - C/m*(x[3]/x[1])^2)
    L = L1 + 2*C / x[1]

    fx = -4*C*C/m/L*x[2]*(x[3]^2)/(x[1]^4) + 2*R*C/m/L*(x[3]^2)/(x[1]^2) + 
        2*C/m*x[2]*(x[3]^2)/(x[1]^3)
    gx = -2*C/m/L*x[3]/(x[1]^2)
    p = 4.6022850*ones(Float32, 3) # p = 4.6022850 just grazes the pedestal with 
    #                              #          voltage_limit = 165.54182708519832
    # p = [1, 5, 10.]
    k = [prod(p), p[1]*p[2]+p[1]*p[3]+p[2]*p[3], sum(p)]
    w = dot(-k, ξ)
    return 1/gx*(w - fx)
end

function curr_cmd_controller(extra_state::AbstractArray, q::AbstractArray, u::AbstractArray)
    m = par[:m]
    g = par[:g]
    C = par[:C]
    L1 = par[:L1]
    R = par[:R]
    qd = par[:qd]

    x = (q[1], u[1], extra_state[1])

    v = -25.0*(x[1]-qd) - 10.0*x[2]
    i_ref = x[1]*sqrt( m/C*(g - v) )
    return R * i_ref - 2 * C * (x[2]*x[3] / x[1] / x[1]) # second term cancels out the nonlinear dynamics
end

function extradynamics(extra_state::AbstractArray, q::AbstractArray, u::AbstractArray)
    m = par[:m]
    g = par[:g]
    C = par[:C]
    L1 = par[:L1]
    R = par[:R]
    qd = par[:qd]

    L = L1 + 2*C / q[1]
    voltage_limit = 165.54182708519832      # This is the minimum voltage needed to pull the ball from δ[2] = 0.2
    
    control_input = clamp(fl_controller(extra_state, q, u), -voltage_limit, voltage_limit)
    # control_input = clamp(curr_cmd_controller(extra_state, q, u), -voltage_limit, voltage_limit)

    return (-R/L * extra_state[1] + 2 * C / L * (u[1]*extra_state[1] / q[1] / q[1]) 
        + 1/L*control_input) * ones(Float32, 1)
end

maglev = Integrator(gap, dynamics, q0, u0, extra_state0, Δt=par[:Δt], ε=par[:ε])
integrate(maglev, extradynamics, par[:tf])


fig = figure(1, clear=true, figsize=(12.80,12.80))
fig.add_subplot(2,2,1)
fig.add_subplot(2,2,2)
fig.add_subplot(2,2,3)
fig.add_subplot(2,2,4)
fig.axes[1].clear()
fig.axes[1].plot(maglev.t, getindex.(maglev.q,1), "b", linewidth=2, label="ball position")
fig.axes[2].plot(maglev.t, getindex.(maglev.u,1), "k", linewidth=2, label="ball velocity")
fig.axes[3].plot(maglev.t, getindex.(maglev.extra_state,1), "m", linewidth=2, label="current")
# fig.axes[4].plot(maglev.t[2:end], getindex.(maglev.Λ,1)./maglev.Δt, "m", linewidth=2, label="contact force top")
# fig.axes[4].plot(maglev.t[2:end], getindex.(maglev.Λ,2)./maglev.Δt, "g:", linewidth=2, label="current force bottom")
fig.axes[4].plot(maglev.t, par[:g] .- par[:C]/par[:m]*(getindex.(maglev.extra_state,1)./getindex.(maglev.q,1)).^2, 
    "g", linewidth=2, label="contact force top")

fig.axes[1].set_ylabel(L"$z(t)$", fontsize=18)
fig.axes[2].set_ylabel(L"$\dot{z}(t)$", fontsize=18)
fig.axes[3].set_ylabel(L"$i(t)$", fontsize=18)
# fig.axes[4].set_ylabel(L"$\Lambda(t)$", fontsize=18)
fig.axes[4].set_ylabel(L"$\xi_3(t)$", fontsize=18)

# fig.axes[4].legend(loc="best", fontsize="x-large")

for i = 1:4
    fig.axes[i].set_xlabel(L"$t$", fontsize=18)
    # fig.axes[i].autoscale_view()
    fig.axes[i].grid()
end
fig.suptitle("Magnetic levitation", fontsize=18)
fig.tight_layout()


# fig = figure(2, clear=true, figsize=(12.80,12.80))
# fig.add_subplot(1,1,1)
# fig.axes[1].clear()
# fig.axes[1].plot(maglev.t, 
#     getindex.(maglev.q,1) .- getindex.(maglev.u,1) .+  getindex.(maglev.extra_state,1),
#     linewidth=2,
#     label="eigvector" )