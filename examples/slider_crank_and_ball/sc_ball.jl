"""
Slider crank with ball
q = (θ1,θ2,x3,x4), u = (θ1dot,θ2dot,x3dot,x4dot),
g(q) = x4 - (x3 + w)
W = del g / del q = [0; 0; -1; 1]
"""

using LinearAlgebra
using NLsolve
using PyPlot
using TimeStepping

const l1 = 0.5
const l2 = 0.75
mass = 1.0
g = 9.81
qA = [π/3; -π/4; 0.0; 1.0]
uA = [0.0; 0.0; 0.0; -1.0]

function loop_fixed(q::Vector)
    return [
        l1*cos(qA[1]) + l2*cos(q[1]) - q[2],
        l1*sin(qA[1]) + l2*sin(q[1])
    ]
end

sol = nlsolve(q->loop_fixed(q[1:2]), [-π/4, 1.0], autodiff=:forward)
qA[2:3] = sol.zero

function loop(q::Vector)
    return [
        l1*cos(q[1]) + l2*cos(q[2]) - q[3],
        l1*sin(q[1]) + l2*sin(q[2])
    ]
end

function jac(q::Vector)
    jac = zeros(Float64, 2, 4)
    jac[1,1] = -l1*sin(q[1])
    jac[1,2] = -l2*sin(q[2])
    jac[1,3] = -1
    jac[2,1] = l1*cos(q[1])
    jac[2,2] = l2*cos(q[2])
    return jac
end


function gap(q::Vector, u::Vector)
    W = zeros(4,1)
    W[3,1] = -1.0
    W[4,1] = 1.0
    return [q[4]-q[3]-0.1], W
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