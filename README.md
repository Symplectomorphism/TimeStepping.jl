# TimeStepping

[![Coverage](https://codecov.io/gh/Symplectomorphism/TimeStepping.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Symplectomorphism/TimeStepping.jl)

## Basic Usage

1. Define the _gap_ function with the signature
> gap(q::Vector, u::Vector)\
> return gap, force_matrix

2. Define the _dynamics_ function with the signature
> dynamics(q::Vector, u::Vector)\
> return M, h

3. Initialize the system by providing the gap and dynamics functions, and the initial conditions.
> bouncing_ball = Integrator(gap, dynamics, qA, uA; Δt=1e-3)

4. Integrate the system until final_time
> integrate(bouncing_ball, 2.0)

### Examples
Examples are provided within the _examples_ folder.
