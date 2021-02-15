module TimeStepping

using ForwardDiff
using JuMP
using PATHSolver
using Ipopt
using LinearAlgebra
using Revise

include("extra_file.jl")
include("moreau.jl")

export my_f, derivative_of_my_f

end
