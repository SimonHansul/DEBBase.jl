module DEBBase

using DEBParamStructs
using DoseResponse

using Parameters
using ComponentArrays
using DifferentialEquations
using DocStringExtensions
using DataFrames
using PrecompileTools
using StaticArrays

# FIXME: "using DEBBase" takes 140 seconds and 2 GB allocs...

include("Structures.jl")
include("IO.jl")
include("ModelFunctions.jl")
include("Simulators.jl")
include("ImpliedTraits.jl")

@compile_workload begin
    sol = simulator(BaseParamCollection())
end

export AbstractParams,
AbstractStatevars,
AbstractParamCollection,
GlobalBaseParams, 
GlobalBaseStatevars, 
DEBBaseParams, 
BaseParamCollection,
simulator,
isolate_pmoas!

end # module DEBBase
