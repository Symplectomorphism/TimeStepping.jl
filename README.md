# TimeStepping

[![Coverage](https://codecov.io/gh/Symplectomorphism/TimeStepping.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Symplectomorphism/TimeStepping.jl)

## Basic Usage

1. Define the _gap_ function with the signature
> gap(q::Vector, u::Vector)\
> return gap, force_matrix

2. Define the _dynamics_ function with the signature
> dynamics(q::Vector, u::Vector)\
> return M, h

3. If there are bilateral constraints, define the holonomic constraint functions and their derivatives (Jacobian and Jacobian_dot).
> loop(q::Vector), jac(q::Vector), jacdot(q::Vector, u::Vector)\
> return ϕ, J, Jdot

4. Initialize the system by providing the gap and dynamics functions, and the initial conditions (for holonomic constraints, provide them, too).
> bouncing_ball = Integrator(gap, dynamics, qA, uA; Δt=1e-3)
> sc_ball = Integrator(gap, dynamics, loop, jac, jacdot, q0, u0; Δt=1e-3)

5. Integrate the system until final_time
> integrate(bouncing_ball, final_time)

### Examples
Examples are provided within the _examples_ folder.
