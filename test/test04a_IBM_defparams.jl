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

# OPTION 1:
# Make an exception for Idot_max_rel and add it as separate field to 


DEBBaseParams()

"""
DEBBase Agent. <br>
Each agent owns a reference to its associated parameter collection.
"""
mutable struct BaseAgent
    p::Base.RefValue{BaseParamCollection} # reference to the paramter collection
    
    function BaseAgent(p::A) where A <: AbstractParamCollection # generate agent from paramter collection
        ref = Ref(p)
        a = new()
        a.p = ref
        a.u = initialize_statevars(a.p)
        a.du = similar(a.u)
        return a
    end
end

"""
    individual_variability!(p::Ref{AbstractParams})
Induce individual to DEB parameters via zoom factor `Z`. 
`Z` is sampled from the corresponding distribution given in `p` and assumed to represent a ratio between maximum structurel *masses* (not lengths), 
so that the surface area-specific ingestion rate `Idot_max_rel` scales with `Z^(1/3)` and parameters which represent masses or energy pools scales with `Z`.
"""
function individual_variability!(p::Ref{AbstractParams})
    Z_a = rand(p.x.Z)
    p.x.Idot_max_rel = p.x.Idot_max_rel_mean * Z_a^(1/3)
    p.x.Idot_max_rel_emb = p.x.Idot_max_rel_emb_mean * Z_a^(1/3)

    for param in p.x.propagate_zoom
        setproperty!()
    end
end


a = BaseAgent(BaseParamCollection())

# Idot_max_rel => Distribution(Idot_max_rel_mean, Idot_max_rel_sd)
# 

