module TimeStepping

using ForwardDiff
using JuMP
using PATHSolver
using Ipopt
using LinearAlgebra

include("extra_file.jl")

export my_f, derivative_of_my_f

end
