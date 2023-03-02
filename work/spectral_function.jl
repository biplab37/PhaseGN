using PhaseGN

p = Parameters(κ=0.01)

sf(ω,q) = structure_factor_ϕ(0.78,0.0,ω,q,p)

using Plots
contourf(0.0:0.02:1.3,0.0:0.02:1.3,sf)