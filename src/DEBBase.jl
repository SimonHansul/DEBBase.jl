module DEBBase

using Parameters
using DifferentialEquations
using DocStringExtensions
using DataFrames
using PrecompileTools

include("Structures.jl")
include("IO.jl")
include("ModelFunctions.jl")
include("Simulators.jl")
include("ImpliedTraits.jl")

export GlobalBaseParams, 
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
