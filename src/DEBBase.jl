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
export relative_response

include("ModelFunctions.jl")

include("Simulators.jl")
export simulator, @replicates

include("ImpliedTraits.jl")
include("ParamHandling.jl")

@compile_workload begin
    # precompile the default simulator
    yhat = simulator(DEBParamCollection())
end


export isolate_pmoas!, sig, clipneg, set_equal!

end # module DEBBase
