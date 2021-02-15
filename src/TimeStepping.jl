module TimeStepping

using DataStructures
using ForwardDiff
using JuMP
using PATHSolver
using Ipopt
using LinearAlgebra
using Revise
import MathOptInterface

include("extra_file.jl")
include("moreau.jl")
include("integrator.jl")

export my_f, derivative_of_my_f
export Moreau, step
export Integrator, integrate

end
