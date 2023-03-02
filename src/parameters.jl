Base.@kwdef mutable struct Parameters
    # parameters for the model.
    Λ::Float64 = 100.0 # cutoff for the regularisation
    κ::Float64 = 0.0
    M::Float64 = 1.0
end

import Base.Broadcast: broadcastable

broadcastable(p::Parameters) = Ref(p)

function Base.show(io::IO,::MIME"text/plain",p::Parameters)
    println(io,"Parameters:")
    println(io,"Λ = ",p.Λ)
    println(io,"κ = ",p.κ)
    println(io,"M = ",p.M)
end

function default_parameters()
    p = Parameters()
end

export Parameters, default_parameters