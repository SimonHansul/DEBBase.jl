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


module ParamStructs
    
    include("ParamStructs/paramstructs.jl") # definition of type hierarchy for parameter structures
    include("ParamStructs/structgeneration.jl") # functions to generate new parameter structures from base params (experimental)

    export AbstractParams, AbstractSpeciesParams, AbstractGlobalParams, AbstractParamCollection
end

module Utils

    using CSV
    using DataFrames
    using OrdinaryDiffEq
    using ComponentArrays
    using StatsBase

    using ..ParamStructs: AbstractParams, AbstractSpeciesParams, AbstractGlobalParams, AbstractParamCollection
    
    include("Utils/utils.jl")
    export skipinf, vectify, which_in, geomrange, diffvec, fround, drop_na, drop_na!, replace_na!, get_treatment_names, lab, read_W3C, ismin, sol_to_df, sol_to_mat

    include("Utils/ioutils.jl")

    include("Utils/inputprocessing.jl")
    export into!, into, isolate_pmoas!, isolate_pmoas, set_equal!

    include("Utils/outputprocessing.jl")
    export relative_response, idcol!

end

module DoseResponse
    include("DoseResponse/doseresponse.jl")
    export LL2, LL2h, LL2M
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

    using ..ParamStructs: AbstractParams, AbstractSpeciesParams, AbstractGlobalParams, AbstractParamCollection
    using ..DoseResponse: LL2, LL2M, LL2h
    using ..Utils: sol_to_df, sol_to_mat

    include("DEBODE/paramstructs.jl")
    export AbstractABM, GlobalParams, GlobalBaseStatevars, SpeciesParams, DEBParamCollection, ODEAgentParams

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
    using ..ParamStructs: AbstractParams

    include("ABC/structs.jl")
    export Priors, SMCResult, ppc

    include("ABC/paramhandling.jl")
    export assign!, getparam

    include("ABC/sampling.jl")
    export rand!, posterior_sample!

    include("ABC/initialization.jl") # initialization of 
    export deftruncnorm, deflognorm

    include("ABC/smc.jl") # parameter estimation using sequential monte carlo approximate bayesian computation
    export SMC
    
    include("ABC/evaluation.jl") # evaluation of smc output
    export summarize_accepted, ppc
      
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
