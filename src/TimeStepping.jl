module TimeStepping

using JuMP
using PATHSolver
using Ipopt

include("extra_file.jl")

export my_f, derivative_of_my_f

end
