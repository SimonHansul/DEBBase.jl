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
export GlobalParams, GlobalBaseStatevars, SpeciesParams, DEBParamCollection

include("IO.jl")
export setproperty!, isolate_pmoas!, set_equal!

export relative_response

include("ModelFunctions.jl")
export sig, clipneg

include("Simulators.jl")
export simulator, @replicates

include("ImpliedTraits.jl")

@compile_workload begin
    # precompile the default simulator
    yhat = simulator(DEBParamCollection())
end

end # module DEBBase
