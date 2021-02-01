mutable struct Dynamics
    M::Array{Float64, 2}
    h::Array{Float64, 1}
    Δt::Float64
    qA::Array{Float64, 1}
    uA::Array{Float64, 1}
    qM::Array{Float64, 1}
    qE::Array{Float64, 1}
    uE::Array{Float64, 1}
    W::Array{Float64, 2}            # Jacobian transpose in the normal direction
    Λ::Array{Float64, 1}            # Normal contact force
    g::Function                     # Gap function
    ε::Float64                      # Coefficient of restitution (normal dir.)
end


function Dynamics(M::Array{Float64, 2}, h::Array{Float64, 1})

end