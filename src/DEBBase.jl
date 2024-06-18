module DEBBase

using Reexport
using Parameters
using ComponentArrays
using OrdinaryDiffEq
using Distributions
using DocStringExtensions
using DataFrames
using PrecompileTools
using StaticArrays
using StatsBase

module ParamStructs
    include("ParamStructs/paramstructs.jl") # definition of type hierarchy for parameter structures
    include("ParamStructs/structgeneration.jl") # functions to generate new parameter structures from base params (experimental)
end

module DoseResponse
    include("DoseResponse/doseresponse.jl")
end

module DEB
    using Parameters
    using ComponentArrays
    using OrdinaryDiffEq
    using Distributions
    using DocStringExtensions
    using DataFrames
    using PrecompileTools
    using StaticArrays
    using StatsBase

    include("structs.jl")
    export AbstractABM, AbstractAgent, GlobalParams, GlobalBaseStatevars, SpeciesParams, DEBParamCollection, AgentParams

    include("IO.jl")
    export setproperty!, isolate_pmoas!, set_equal!

    export relative_response

    include("derivatives.jl")
    export sig, clipneg

    include("simulators.jl")
    export init_substates_agent, init_substates_global, abstractsimulator, returntypes, simulator, @replicates

    include("traits.jl")

    @compile_workload begin
        # precompile the default simulator
        yhat = simulator(DEBParamCollection())
    end
end

module ABC
    using DataFrames
    using SHUtils
    using DocStringExtensions
    using StatsBase
    using KernelDensity
    using Distributions
    using Parameters
    using RecipesBase
    import Base:rand
    using Base.Threads
    using Dates
    using Random

    include("ABC/structs.jl")
    export Priors, get, SMCResult, opc

    include("ABC/paramhandling.jl")
    export assign!, getparam

    include("ABC/sampling.jl")
    export rand!, posterior_sample!

    include("ABC/initialization.jl")
    export deftruncnorm, deflognorm

    include("ABC/evaluation.jl")
    export summarize_accpeted, ppc

    include("ABC/smc.jl")
    export SMC

    include("ABC/recipes.jl")
end

end # module DEBBase
