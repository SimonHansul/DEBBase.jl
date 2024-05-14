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
export AbstractABM, AbstractAgent, GlobalParams, GlobalBaseStatevars, SpeciesParams, DEBParamCollection, AgentParams

include("IO.jl")
export setproperty!, isolate_pmoas!, set_equal!

export relative_response

include("ModelFunctions.jl")
export sig, clipneg

include("Simulators.jl")
export init_substates_agent, init_substates_global, abstractsimulator, returntypes, simulator, @replicates

include("ImpliedTraits.jl")

@compile_workload begin
    # precompile the default simulator
    yhat = simulator(DEBParamCollection())
end

end # module DEBBase
