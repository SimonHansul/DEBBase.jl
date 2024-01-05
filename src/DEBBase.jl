module DEBBase

using Parameters
using ComponentArrays
using DoseResponse
using DifferentialEquations
using DocStringExtensions
using DataFrames
using PrecompileTools

# FIXME: "using DEBBase" takes 140 seconds and 2 GB allocs...

include("Structures.jl")
include("IO.jl")
include("ModelFunctions.jl")
include("Simulators.jl")
include("ImpliedTraits.jl")

export AbstractParams,
AbstractStatevars,
GlobalBaseParams, 
GlobalBaseStatevars, 
DEBBaseParams, 
BaseParamCollection,
simulator,
isolate_pmoas!

@compile_workload begin
    sol = simulator(BaseParams())
end

end # module DEBBase
