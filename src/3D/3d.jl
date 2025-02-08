module ThreeD

using Cubature, PolyLog
# using UsefulFunctions: PVintegral
using NLsolve: nlsolve

include("../parameters.jl")
include("../helperfunction.jl")

include("custom_types.jl")
include("gapeqns3d.jl")

end
