module DEBBase

using Reexport
@reexport using DEBParamStructs
using DoseResponse
using Parameters
using ComponentArrays
using OrdinaryDiffEq
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
    # precompile the default simulator
    yhat = simulator(DEBParamCollection())
end

export GlobalParams, 
GlobalBaseStatevars, 
SpeciesParams, 
DEBParamCollection,
isolate_pmoas!,
sig,
clipneg,
relative_response,
set_equal!,
simulator,
@replicates,
@sweep

end # module DEBBase
