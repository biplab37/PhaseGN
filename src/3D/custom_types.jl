Base.@kwdef mutable struct Couplings3D
    Gs::Float64 = 
    Gv::Float64
    Gd::Float64
end

Base.@kwdef mutable struct Parameters3D
    cutoff::Float64 = 0.64
    coupling::Couplings3D = Couplings3D()
    m0::Float64 = 0.0054
end