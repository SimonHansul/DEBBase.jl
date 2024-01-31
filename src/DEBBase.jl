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

include("Structures.jl")
include("IO.jl")
include("ModelFunctions.jl")
include("Simulators.jl")
include("ImpliedTraits.jl")

@compile_workload begin
    sol = simulator(BaseParamCollection())
end

export GlobalBaseParams, 
GlobalBaseStatevars, 
DEBBaseParams, 
BaseParamCollection,
simulator,
isolate_pmoas!,
sig,
clipneg

end # module DEBBase
