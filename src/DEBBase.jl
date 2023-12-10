module DEBBase

using Parameters
using DifferentialEquations
using DocStringExtensions

include("IO.jl")
include("ModelFunctions.jl")
include("Simulators.jl")
include("ImpliedTraits.jl")

export GlobalParams, DEBParams, DEB!

end # module DEBBase
