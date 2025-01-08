Base.@kwdef mutable struct Couplings3D
    Gs::Float64 = 2.134 / 0.639^2
    Gv::Float64 = Gs / 2
    Gd::Float64 = 0.75 * Gs
end

Base.@kwdef mutable struct Parameters3D
    cutoff::Float64 = 0.639
    m0::Float64 = 0.0055
    coupling::Couplings3D = Couplings3D()
end

broadcastable(p::Parameters3D) = Ref(p)
broadcastable(g::Couplings3D) = Ref(g)

function Base.show(io::IO, ::MIME"text/plain", p::Parameters3D)
    println(io, "Parameters:")
    println(io, "Λ = ", p.cutoff)
    println(io, "m0 = ", p.m0)
    println(io, "Gs = ", p.coupling.Gs)
    println(io, "Gv = ", p.coupling.Gv)
    println(io, "Gd = ", p.coupling.Gd)
end

function Base.show(io::IO, ::MIME"text/plain", p::Couplings3D)
    println(io, "Couplings:")
    println(io, "Gs = ", p.Gs)
    println(io, "Gv = ", p.Gv)
    println(io, "Gd = ", p.Gd)
end
