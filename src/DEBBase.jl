module DEBBase

using Parameters
using ComponentArrays
using DoseResponse
using DifferentialEquations
using DocStringExtensions
using DataFrames
using PrecompileTools

const G::Int64 = 1
const M::Int64 = 2
const A::Int64 = 3
const R::Int64 = 4
const H::Int64 = 5

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
DEBBaseOrganism,
simulator

@compile_workload begin
    glb = GlobalBaseParams()
    anm = DEBBaseParams()
    sol = simulator(glb, anm)
end

end # module DEBBase
