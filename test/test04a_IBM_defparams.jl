@info("Loading packages")
@time begin
    using Pkg; Pkg.activate("test")
    using Test
    using Plots, StatsPlots, Plots.Measures
    using SHUtils
    using DataFrames, DataFramesMeta
    using ProgressMeter
    default(leg = false, lw = 1.5)
    TAG = splitpath(@__FILE__)[end] |> x -> split(x, ".")[1] |> String    

    using Revise
    using DEBBase
end
using Parameters, DEBParamStructs, DEBBase

#=
How to implement hierarchical models in an IBM context?

We have

    - parameters which are common across individuals of a species, `pcmn::Ref{AbstractParamCollection}`, given as a reference to the param struct
    - parameters which are specific to an individual, `pown::ComponentVector`, newly assigned for each individual

The model functions take `p` as argument. 
How do we incorporate `pown` then? 

**Solution**: 
- Always assume hierarchical model. 
    ~~- All model functions take `pcmn` and `pown` as arguments. 
    ~~- The default Zoom factor is Dirac(1.).
    - simulator() takes care of initializing `pown`
=#

DEBBaseParams()

"""
DEBBase Agent. <br>
Each agent owns a reference to its associated parameter collection.
"""
mutable struct BaseAgent
    p::Base.RefValue{BaseParamCollection} # reference to the paramter collection
    u::CompositeVector
    du::CompositeVector

    function BaseAgent(p::A) where A <: AbstractParamCollection # generate agent from paramter collection
        ref = Ref(p)
        a = new()
        a.p = ref
        a.u = initialize_statevars(a.p)
        a.du = similar(a.u)
        return a
    end
end


a = BaseAgent(BaseParamCollection())

# Idot_max_rel => Distribution(Idot_max_rel_mean, Idot_max_rel_sd)
# 

