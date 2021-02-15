# TimeStepping

[![Coverage](https://codecov.io/gh/Symplectomorphism/TimeStepping.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Symplectomorphism/TimeStepping.jl)

## Example: Bouncing ball

"""
Bouncing ball
q = (x,z), u = (xdot, zdot),
g(q) = q[2]
W = del g / del q = [0; 1]
M = [m 0; 0 m]
h = [0 m*g]
"""

using LinearAlgebra
using PyPlot
using TimeStepping

mass = 1.0
g = 9.81
qA = [0.0; 1.0]
uA = [0.0; -0.5]

function gap(q::Vector, u::Vector)
    W = zeros(2,1)
    W[1,1] = 0.0
    W[2,1] = 1.0
    return [q[2]], W
end

function dynamics(q::Vector, u::Vector)
    return diagm(ones(2)*mass), [0; -mass*g]
end

bouncing_ball = Integrator(gap, dynamics, qA, uA; Δt=1e-3)
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
