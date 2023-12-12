module DEBBase

using Parameters
using DifferentialEquations
using DocStringExtensions

include("Structures.jl")
include("IO.jl")
include("ModelFunctions.jl")
include("Simulators.jl")
include("ImpliedTraits.jl")

export GlobalBaseParams, 
GlobalBaseStatevars, 
DEBBaseParams, 
DEBBaseStatevars,
DEBBaseOrganism

end # module DEBBase
