module DEBBase

using Reexport
using Parameters
using ComponentArrays
using OrdinaryDiffEq
using Distributions
using DataFrames
using PrecompileTools
using StaticArrays
using StatsBase
using NamedTupleTools
using Base.Threads
using Random
using YAML
using DataStructures

 # definition of type hierarchy for parameter structures 
module ParamStructs
    include("ParamStructs/paramstructs.jl")
    export AbstractParams, AbstractSpeciesParams, AbstractGlobalParams, AbstractParamCollection
end

module DoseResponse
    include("DoseResponse/doseresponse.jl")
    export LL2, LL2h, LL2M
end

module Utils

    using CSV
    using DataFrames
    using OrdinaryDiffEq
    using ComponentArrays
    using StatsBase
    using YAML
    using NamedTupleTools
    using DataStructures
    using Distributions

    using ..DoseResponse

    # functions from sister modules have to be imported explictly
    using ..ParamStructs: AbstractParams, AbstractSpeciesParams, AbstractGlobalParams, AbstractParamCollection

    include("Utils/utils.jl")
    export skipinf, vectify, which_in, geomrange, diffvec, fround, drop_na, drop_na!, replace_na!, get_treatment_names, lab, read_W3C, ismin, sol_to_df, sol_to_mat

    include("Utils/ioutils.jl")
    export load_config, save_config

    include("Utils/inputprocessing.jl")
    export into!, into, isolate_pmoas!, isolate_pmoas, set_equal!, getstat, C2K

    include("Utils/outputprocessing.jl")
    export relative_response, idcol!
end



module DEBODE

    using Parameters
    using ComponentArrays
    using OrdinaryDiffEq
    using Distributions
    using DataFrames
    using PrecompileTools
    using StaticArrays
    using StatsBase
    using Base.Threads

    using ..ParamStructs: AbstractParams, AbstractSpeciesParams, AbstractGlobalParams, AbstractParamCollection
    using ..DoseResponse: LL2, LL2M, LL2h
    using ..Utils: sol_to_df, sol_to_mat

    include("DEBODE/paramstructs.jl")
    export AbstractABM, GlobalParams, GlobalBaseStatevars, SpeciesParams, Params, AgentParams

    include("DEBODE/derivatives.jl")
    export sig, clipneg

    include("DEBODE/statevars.jl")
    export initialize_statevars

    include("DEBODE/events.jl")

    include("DEBODE/simulators.jl")
    export @replicates, replicates, treplicates, exposure

    include("DEBODE/traits.jl")

    @compile_workload begin
        # precompile the default simulator
        for _ in 1:10
            sim = simulator(Params())
        end
    end
end

module DEBABM

    using Parameters
    using ComponentArrays
    using NamedTupleTools
    using DataFrames
    using Base.Threads
    using Random


    using ..DEBODE: AgentParams, GlobalParams, SpeciesParams, Params # import ODE params
    using ..DEBODE: DEBODE_global!, DEBODE_agent_IA! # import ODE derivatives
    using ..DEBODE: initialize_agent_statevars, initialize_global_statevars # import ODE statevars
    using ..DEBODE: condition_juvenile, condition_adult, effect_juvenile!, effect_adult! # import ODE lifestage definitions
    using ..DEBODE: sig
    using ..DoseResponse: LL2h, LL2
    using ..ParamStructs: AbstractParamCollection, AbstractParams # import paramstructs

    include("DEBABM/paramstructs.jl")
    export ABMParamCollection
    include("DEBABM/agents.jl")
    include("DEBABM/model.jl")
    include("DEBABM/schedules.jl")
    include("DEBABM/simulators.jl")

    export AbstractDEBAgent, AbstractDEBABM, DEBAgent, DEBABM, agent_record_to_df
end

"""
Submodule for parameter estimation using approximate bayesian computation.
"""
module ABC

    using DataFrames
    using StatsBase
    using KernelDensity
    using Distributions
    using Parameters
    using RecipesBase
    import Base:rand
    using Base.Threads
    using Dates
    using Random

    using ..Utils: fround
    using ..ParamStructs: AbstractParams, AbstractParamCollection

    include("ABC/priors.jl")
    export Priors, ppc

    include("ABC/paramhandling.jl")
    export assign!, getparam

    include("ABC/smc.jl") # parameter estimation using sequential monte carlo approximate bayesian computation
    export SMC, SMCResult

    include("ABC/sampling.jl")
    export rand!, posterior_sample!

    include("ABC/initialization.jl") # initialization of 
    export deftruncnorm, deflognorm

    include("ABC/evaluation.jl") # evaluation of smc output
    export summarize_accepted, ppc

    include("ABC/imanconover.jl") # implementation of Iman-Conover algorithm
      
    include("ABC/recipes.jl")
end

"""
Recipes and convenience functions for plotting model output.
"""
module Figures

    using StatsBase
    using Reexport
    using RecipesBase

    include("Figures/figutils.jl")
    export  gridylabel, gridxlabel, gridyattr, gridxattr

    include("Figures/recipes.jl")

    export rugplot, lineplot, groupedlineplot
end



end # module DEBBase
