module ThreeD

using Cubature, PolyLog
# using UsefulFunctions: PVintegral
using NLsolve: nlsolve

include("../parameters.jl")
include("../helperfunction.jl")

include("gapeqns3d.jl")
include("custom_types.jl")

end
