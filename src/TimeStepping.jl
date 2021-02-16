module TimeStepping

using DataStructures
using ForwardDiff
using JuMP
using PATHSolver
using Mosek
using MosekTools
using Ipopt
using LinearAlgebra
using Revise
import MathOptInterface

include("extra_file.jl")
include("moreau.jl")
include("integrator.jl")

export my_f, derivative_of_my_f
export Moreau, step, step_constrained
export Integrator, integrate

end
