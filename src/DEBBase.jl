module DEBBase

using Reexport
@reexport using DEBParamStructs
@reexport using DoseResponse
using ComponentArrays, StaticArrays
using OrdinaryDiffEq, Agents
using Distributions, StatsBase, Random
using DataFrames
using PrecompileTools

include("Solvers.jl")

include("Structures.jl")
export AbstractABM, AbstractSpeciesParams, ABM, BaseAgent, GlobalParams, GlobalBaseStatevars, SpeciesParams, DEBParamCollection, AgentParams

include("Initialize.jl")
export init_substates_agent, init_substates_global, initialize_statevars, initialize_statevars!, initialize_agents!

include("IO.jl")
export setproperty!, isolate_pmoas!, set_equal!, relative_response

include("ModelFunctions.jl")
export sig, clipneg

include("Simulators.jl")
export init_substates_agent, init_substates_global, abstractsimulator, returntypes, simulator, @replicates

include("ImpliedTraits.jl")

@compile_workload begin
    # precompile the default simulator
    #yhat = simulator(DEBParamCollection())
end

end # module DEBBase
