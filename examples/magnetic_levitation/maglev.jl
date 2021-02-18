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

const m = 0.01187
const g = 9.81
const C = 1.4e-4
const L1 = 0.65
const R = 28.7
final_time = 2.0
q0 = convert.(Float32, [0.19])           # Start the ball on the pedestal,
u0 = convert.(Float32, [0.0])           # with zero velocity
extra_state0 = convert.(Float32, [0.0]) # and zero current.

function gap(q::AbstractArray, u::AbstractArray)
    return [q[1].-0.05, 0.20.-q[1]], convert.(Float32, [1.0 -1])
end

function dynamics(q::AbstractArray, u::AbstractArray, extra_state::AbstractArray)
    return convert.(Float32, m)*Matrix(1.0f0*I, 1, 1), convert.(Float32, m*g .- C*(extra_state[1]./q[1]).^2)*ones(Float32, 1)
end

function extradynamics(extra_state::AbstractArray, q::AbstractArray, u::AbstractArray)
    L = L1 + 2*C / q[1]
    ## Formulate the controller (could be put in another function)
    voltage_limit = 1.0
    v = -25.0*(q[1]-0.1) - 10.0*u[1]
    v = clamp(v, -voltage_limit, voltage_limit)
    i_ref = q[1]*sqrt( m/C*(g - v) )
    control_input = R * i_ref - 2 * C * (u[1]*extra_state[1] / q[1] / q[1]) # second term cancels out the nonlinear dynamics
    ## end controller formulation
    return (-R/L * extra_state[1] + 2 * C / L * (u[1]*extra_state[1] / q[1] / q[1]) + 1/L*control_input) * ones(Float32, 1)
end

maglev = Integrator(gap, dynamics, convert.(Float32, q0), convert.(Float32, u0), 
    convert.(Float32, extra_state0); Δt=convert(Float32, 1e-3))
integrate(maglev, extradynamics, final_time)


fig = figure(1, clear=true, figsize=(12.80,13.85))
fig.add_subplot(2,2,1)
fig.add_subplot(2,2,2)
fig.add_subplot(2,2,3)
fig.add_subplot(2,2,4)
fig.axes[1].clear()
fig.axes[1].plot(maglev.t, getindex.(maglev.q,1), "b", linewidth=2, label="ball position")
fig.axes[2].plot(maglev.t, getindex.(maglev.u,1), "k", linewidth=2, label="ball velocity")
fig.axes[3].plot(maglev.t, getindex.(maglev.extra_state,1), "m", linewidth=2, label="current")
fig.axes[4].plot(maglev.t[2:end], getindex.(maglev.Λ,1), "m", linewidth=2, label="contact force top")
fig.axes[4].plot(maglev.t[2:end], getindex.(maglev.Λ,2), "g:", linewidth=2, label="current force bottom")


fig.axes[1].set_xlabel(L"$t$", fontsize=18)
fig.axes[1].set_ylabel(L"$z(t)$", fontsize=18)
fig.axes[2].set_xlabel(L"$t$", fontsize=18)
fig.axes[2].set_ylabel(L"$\dot{z}(t)$", fontsize=18)
fig.axes[3].set_xlabel(L"$t$", fontsize=18)
fig.axes[3].set_ylabel(L"$i(t)$", fontsize=18)
fig.axes[4].set_xlabel(L"$t$", fontsize=18)
fig.axes[4].set_ylabel(L"$\Lambda(t)$", fontsize=18)

fig.axes[1].set_title("Magnetic levitation", fontsize=18)
fig.axes[2].set_title("Magnetic levitation", fontsize=18)
fig.axes[4].legend(loc="best", fontsize="x-large")

fig.axes[1].autoscale_view()
fig.axes[2].autoscale_view()
fig.axes[3].autoscale_view()
fig.axes[4].autoscale_view()
fig.tight_layout()
