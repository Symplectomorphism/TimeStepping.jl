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
const m1 = 0.25
const m2 = 0.4
const m3 = 0.5
const m4 = 1.0
const J1 = 0.05
const J2 = 0.07
const g = 9.81
q0 = [π/3; -π/4; 0.0; 1.0]
u0 = [0.0; 0.0; 0.0; -0.5]

function loop_fixed(q::Vector)
    return [
        l1*cos(q0[1]) + l2*cos(q[1]) - q[2],
        l1*sin(q0[1]) + l2*sin(q[1])
    ]
end

sol = nlsolve(q->loop_fixed(q[1:2]), [-π/4, 1.0], autodiff=:forward)
q0[2:3] = sol.zero

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

function jacdot(q::Vector, u::Vector)
    jacdot = zeros(Float64, 2, 4)
    jacdot[1,1] = -l1*cos(q[1])*u[1]
    jacdot[1,2] = -l2*cos(q[2])*u[2]
    jacdot[2,1] = -l1*sin(q[1])*u[1]
    jacdot[2,2] = -l2*sin(q[2])*u[2]
    return jacdot
end


function gap(q::Vector, u::Vector)
    W = zeros(4,1)
    W[3,1] = -1.0
    W[4,1] = 1.0
    return [q[4]-q[3]-0.1], W
end

function dynamics(q::Vector, u::Vector)
    M = zeros(Float64,4, 4)
    M[1,1] = 1/4*m1*l1*l1 + m2*l1*l1 + J1
    M[1,2] = M[2,1] = 1/2*m2*l1*l2*cos(q[1]-q[2])
    M[2,2] = 1/4*m2*l2*l2 + J2
    M[3,3] = m3
    M[4,4] = m4

    h = zeros(Float64, 4)
    h[1] = -m2*sin(q[1]-q[2])u[1]*u[2] + 1/2*m2*sin(q[1]-q[2])*u[2]*u[2] + (1/2*m1+m2)*g*l1*cos(q[1])
    h[2] = m2*sin(q[1]-q[2])u[1]*u[2] - 1/2*m2*sin(q[1]-q[2])*u[1]*u[1] + m2*g*l1*cos(q[2])
    h[3] = 0.0
    h[4] = 0.0
    return M, -h
end

# m = Moreau(gap, dynamics, loop, jac, jacdot, q0, u0)
sc_ball = Integrator(gap, dynamics, loop, jac, jacdot, q0, u0; Δt=1e-3)
integrate(sc_ball, 2.0)


fig = figure(1, clear=true)
fig.add_subplot(1,1,1)
fig.axes[1].clear()
fig.axes[1].plot(sc_ball.t, getindex.(sc_ball.q,4), "k", linewidth=2)

fig.axes[1].set_title("Ball motion", fontsize=18)
fig.axes[1].set_xlabel(L"$t$", fontsize=18)
fig.axes[1].set_ylabel(L"$x_4(t)$", fontsize=18)

fig.axes[1].autoscale_view()
fig.tight_layout()