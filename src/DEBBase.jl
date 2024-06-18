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


module Utils
        
    using CSV
    using DataFrames

    include("Utils/utils.jl")

    export skipinf, vectify, which_in, geomrange, diffvec, fround, drop_na, drop_na!, replace_na!, get_treatment_names, lab, read_W3C, ismin
end

module ParamStructs
    include("ParamStructs/paramstructs.jl") # definition of type hierarchy for parameter structures
    include("ParamStructs/structgeneration.jl") # functions to generate new parameter structures from base params (experimental)

    export AbstractParams, AbstractSpeciesParams, AbstractGlobalParams, AbstractParamCollection
end

module DoseResponse
    include("DoseResponse/doseresponse.jl")
end

module DEBODE
    using Parameters
    using ComponentArrays
    using OrdinaryDiffEq
    using Distributions
    using DocStringExtensions
    using DataFrames
    using PrecompileTools
    using StaticArrays
    using StatsBase

    using ..ParamStructs
    using ..DoseResponse

    include("DEBODE/paramstructs.jl")
    export AbstractABM, GlobalParams, GlobalBaseStatevars, SpeciesParams, DEBParamCollection, AgentParams

    include("DEBODE/IO.jl")
    export setproperty!, isolate_pmoas!, set_equal!, relative_response

    include("DEBODE/derivatives.jl")
    export sig, clipneg

    include("DEBODE/statevars.jl")

    include("DEBODE/simulators.jl")
    export initialize_statevars, abstractsimulator, simulator, @replicates

    include("DEBODE/traits.jl")

    @compile_workload begin
        # precompile the default simulator
        yhat = simulator(DEBParamCollection())
    end
end

module ABC
    using DataFrames
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

    using ..Utils

    include("ABC/structs.jl")
    export Priors, get, SMCResult, opc

    include("ABC/paramhandling.jl")
    export assign!, getparam

    include("ABC/sampling.jl")
    export rand!, posterior_sample!

    include("ABC/initialization.jl")
    export deftruncnorm, deflognorm

    include("ABC/evaluation.jl")
    export summarize_accepted, ppc

    include("ABC/smc.jl")
    export SMC

    include("ABC/recipes.jl")
end

end # module DEBBase
