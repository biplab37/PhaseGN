Base.@kwdef mutable struct Parameters
    # parameters for the model.
    Λ::Float64 = 100.0 # cutoff for the regularisation
    κ::Float64 = 0.01 # Chiral limit κ->0 has some problem with the regularisation.
    M::Float64 = 1.0
end

# import Base.Broadcast: broadcastable
#
# broadcastable(p::Parameters) = Ref(p)

function Base.show(io::IO, ::MIME"text/plain", p::Parameters)
    println(io, "Parameters:")
    println(io, "Λ = ", p.Λ)
    println(io, "κ = ", p.κ)
    println(io, "M = ", p.M)
end

default_parameters() = Parameters()

export Parameters, default_parameters
