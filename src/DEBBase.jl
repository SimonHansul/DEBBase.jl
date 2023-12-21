module DEBBase

using Parameters
using DoseResponse
using DifferentialEquations
using DocStringExtensions
using DataFrames
using PrecompileTools

include("Structures.jl")
include("ModelFunctions.jl")
include("Simulators.jl")
include("ImpliedTraits.jl")

export AbstractParams,
AbstractStatevars,
GlobalBaseParams, 
GlobalBaseStatevars, 
DEBBaseParams, 
DEBBaseStatevars,
DEBBaseOrganism

@compile_workload begin
    glb = GlobalBaseParams()
    anm = DEBBaseParams()
    sol = DEBBase.run_model(glb, anm)
end

end # module DEBBase
