"""
    PhaseGN

This module contains many of the functions related to the calculation of various quantities in 2+1 dimensional GN model.
"""
module PhaseGN

using Cubature

include("parameters.jl")
include("helperfunction.jl")
include("masses.jl")
include("meanfield.jl")
include("sigma.jl")
include("phi.jl")
include("pressure.jl")
include("width.jl")
include("ext_momentum.jl")

end