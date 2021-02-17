"""
Three bouncing balls
q = (z1,z2,z3), u = (z1dot, z2dot, z3dot),
g(q) = (q[1],q[2],q[3])
W = del g / del q = I (3 x 3 identity)
M = [m1 0 0; 0 m2 0; 0 0 m3]
h = [-m1*g; -m2*g; -m3*g]
"""

using LinearAlgebra
using PyPlot
using TimeStepping

const m1 = m2 = m3 = 1.0
const g = 9.81
final_time = 2.0
q0 = [0.8; 1.0; 1.2]
u0 = [0.5; 2.0; -0.2]

function gap(q::AbstractArray, u::AbstractArray)
    return q, Matrix(1.0f0*I, 3, 3)
end

function dynamics(q::AbstractArray, u::AbstractArray)
    return convert.(Float32, [m1 0 0; 0 m2 0; 0 0 m3]), convert.(Float32, [-m1*g; -m2*g; -m3*g])
end

three_bb = Integrator(gap, dynamics, convert.(Float32, q0), convert.(Float32, u0); Δt=convert(Float32, 1e-3))
integrate(three_bb, final_time)


fig = figure(1, clear=true)
fig.add_subplot(1,1,1)
fig.axes[1].clear()
fig.axes[1].plot(three_bb.t, getindex.(three_bb.q,1), "b", linewidth=2, label="ball 1")
fig.axes[1].plot(three_bb.t, getindex.(three_bb.q,2), "k", linewidth=2, label="ball 2")
fig.axes[1].plot(three_bb.t, getindex.(three_bb.q,3), "m", linewidth=2, label="ball 3")

fig.axes[1].set_title("Three bouncing balls", fontsize=18)
fig.axes[1].set_xlabel(L"$t$", fontsize=18)
fig.axes[1].set_ylabel(L"$z(t)$", fontsize=18)
fig.axes[1].legend(loc="best", fontsize="x-large")

fig.axes[1].autoscale_view()
fig.tight_layout()