module DEBBase

using DEBParamStructs
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

    # precompilation of the @compose workflow
    functions = [ # define ODE system as list of derivative functions
        DEBBase.Idot!,
        DEBBase.Adot!,
        DEBBase.Mdot!,
        DEBBase.Jdot!,
        DEBBase.Sdot!,
        DEBBase.Hdot!,
        DEBBase.H_bdot!,
        DEBBase.Rdot!,
        DEBBase.X_pdot!,
        DEBBase.X_embdot!,
        DEBBase.Ddot!,
        DEBBase.C_Wdot!
    ]

    system! = DEBBase.@compose functions # use @compose to put the functions together
    
    theta = DEBParamCollection()
    theta.glb.t_max = 56.

    yhat = simulator(theta; system = system!)
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
@compose

end # module DEBBase
