module DEBBase

using Parameters
using ComponentArrays
 # TODO: add DoseResponse as dependency again...
 # FIXME: why does DoseResponse error when adding DEBBase as dependency of Pathogens?
#using DoseResponse
using DifferentialEquations
using DocStringExtensions
using DataFrames
using PrecompileTools

# FIXME: "using DEBBase" takes 140 seconds and 2 GB allocs...

include("DRFunctions.jl")
include("Structures.jl")
include("IO.jl")
include("ModelFunctions.jl")
include("Simulators.jl")
include("ImpliedTraits.jl")

export AbstractParams,
AbstractStatevars,
AbstractParamCollection,
GlobalBaseParams, 
GlobalBaseStatevars, 
DEBBaseParams, 
BaseParamCollection,
simulator,
isolate_pmoas!

@compile_workload begin
    sol = simulator(BaseParamCollection())
end

end # module DEBBase
