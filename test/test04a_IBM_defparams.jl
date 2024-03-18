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

# TODO: how to deal with hyperparameters?
# Idot_max_rel is fixed in BaseAgent.p
# But the realized Idot_max_rel is a sample from N(μ_I, σ_I).
# Within du(du, u, p, t), we have to access Idot_max_rel as p.x.deb.Idot_max_rel

# OPTION 1:
# Make an exception for Idot_max_rel and add it as separate field to 

mutable struct BaseHyParams

end

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
        return a
    end
end


a = BaseAgent(BaseParamCollection())

# Idot_max_rel => Distribution(Idot_max_rel_mean, Idot_max_rel_sd)
# 

