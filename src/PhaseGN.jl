"""
    PhaseGN

This module contains many of the functions related to the calculation of various quantities in 2+1 dimensional GN model.
"""
module PhaseGN

using UsefulFunctions: PVintegral, numberB
using Cubature, PolyLog
using NLsolve: nlsolve
using RecipesBase
using Interpolations
using AdaptiveSparseGrids

include("parameters.jl")
include("helperfunction.jl")
include("masses.jl")
include("meanfield.jl")
include("sigma.jl")
include("phi.jl")
include("pressure.jl")
include("width.jl")
include("ext_momentum.jl")
include("external_momentum.jl")
include("phase.jl")
include("structure_factor.jl")
include("plasmon.jl")
include("pressure_fluc.jl")
include("phase_defn.jl")
include("momentum_effect_on_pressure.jl")
include("mott.jl")
include("phases_refactored.jl")
# include("phase_transition_order.jl")

include("3D/3d.jl")
include("Spectral/spectral.jl")
include("Spectral/plasmon.jl")
# include("gridded_calculations.jl")

include("interpolated_functions.jl")
include("Interpolations/interpolations.jl")

include("entropy.jl")

end
