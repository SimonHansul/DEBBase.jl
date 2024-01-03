module DEBBase

using Parameters
using ComponentArrays
using DoseResponse
using DifferentialEquations
using DocStringExtensions
using DataFrames
using PrecompileTools

# FIXME: "using DEBBase" takes 140 seconds and 2 GB allocs...

include("Structures.jl")
include("IO.jl")
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
simulator,
isolate_pmoas!

@compile_workload begin
    glb = GlobalBaseParams()
    anm = DEBBaseParams()
    sol = simulator(glb, anm)
end

end # module DEBBase
