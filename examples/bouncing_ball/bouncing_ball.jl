"""
Bouncing ball
q = (x,z), u = (xdot, zdot),
g(q) = q[2]
W = del g / del q = [0; 1]
M = [m 0; 0 m]
h = [0 -m*g]
"""

using LinearAlgebra
using PyPlot
using TimeStepping

const mass = 1.0
const g = 9.81
qA = [0.0; 1.0]
uA = [1.0; 2.0]

function gap(q::AbstractArray, u::AbstractArray)
    W = zeros(Float32, 2,1)
    W[1,1] = 0.0
    W[2,1] = 1.0
    return [q[2]], W
end

function dynamics(q::AbstractArray, u::AbstractArray)
    return diagm(ones(Float32, 2)*mass), [0; -mass*g]
end

bouncing_ball = Integrator(gap, dynamics, convert.(Float32, qA), 
    convert.(Float32, uA); Δt=1e-3, ε=0.5)
integrate(bouncing_ball, 2.0)


fig = figure(1, clear=true)
fig.add_subplot(1,1,1)
fig.axes[1].clear()
fig.axes[1].plot(bouncing_ball.t, getindex.(bouncing_ball.q,2), "k", linewidth=2)

fig.axes[1].set_title("Bouncing Ball (z-axis)", fontsize=18)
fig.axes[1].set_xlabel(L"$t$", fontsize=18)
fig.axes[1].set_ylabel(L"$z(t)$", fontsize=18)

fig.axes[1].autoscale_view()
fig.tight_layout()