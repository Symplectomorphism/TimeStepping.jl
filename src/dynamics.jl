mutable struct Dynamics
    M::Array{Float64, 2}
    h::Array{Float64, 1}
    Δt::Float64
    qA::Array{Float64, 1}
    uA::Array{Float64, 1}
    qM::Array{Float64, 1}
    qE::Array{Float64, 1}
    uE::Array{Float64, 1}
    W::Array{Float64, 2}
    Λ::Array{Float64, 1}
    g::Function
end


function Dynamics()
    
end