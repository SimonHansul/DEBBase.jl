module DEBBase

using Parameters
using DifferentialEquations
using DocStringExtensions

include("IO.jl")
include("ModelFunctions.jl")
include("Simulators.jl")

export DEB!

end # module DEBBase
