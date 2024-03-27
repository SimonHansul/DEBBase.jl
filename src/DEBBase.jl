module DEBBase

using Reexport
@reexport using DEBParamStructs
using DoseResponse
using Parameters
using ComponentArrays
using DifferentialEquations
using Distributions
using DocStringExtensions
using DataFrames
using PrecompileTools
using StaticArrays
using StatsBase

include("Structures.jl")
include("IO.jl")
include("ModelFunctions.jl")
include("Simulators.jl")
include("ImpliedTraits.jl")
include("ParamHandling.jl")

@compile_workload begin
    sol = simulator(BaseParamCollection())
end

export GlobalBaseParams, 
GlobalBaseStatevars, 
DEBBaseParams, 
BaseParamCollection,
isolate_pmoas!,
sig,
clipneg,
relative_response,
set_equal!,
@replicates

end # module DEBBase
