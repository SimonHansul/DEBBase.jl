module DEBBase

using Parameters
using ComponentArrays
using DoseResponse
using DifferentialEquations
using DocStringExtensions
using DataFrames
using PrecompileTools

include("IO.jl")
include("Structures.jl")
include("ModelFunctions.jl")
include("Simulators.jl")
include("ImpliedTraits.jl")

glb = GlobalBaseParams()
anm = DEBBaseParams()
sol = simulator(glb, anm)

export AbstractParams,
AbstractStatevars,
GlobalBaseParams, 
GlobalBaseStatevars, 
DEBBaseParams, 
DEBBaseStatevars,
DEBBaseOrganism,
simulator

@compile_workload begin
    glb = GlobalBaseParams()
    anm = DEBBaseParams()
    sol = simulator(glb, anm)
end

end # module DEBBase
