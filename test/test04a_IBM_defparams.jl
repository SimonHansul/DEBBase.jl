@time "Loading packages" begin
    using Pkg; Pkg.activate("test")
    using Test
    using Plots, StatsPlots, Plots.Measures
    using SHUtils
    using DataFrames, DataFramesMeta
    using ProgressMeter
    using Distributions
    default(leg = false, lw = 1.5)
    TAG = splitpath(@__FILE__)[end] |> x -> split(x, ".")[1] |> String    

    using Revise
    using DEBBase
    using Parameters, DEBParamStructs
end

begin
    plt = plot(layout = (1,2))

    y = DEBBase.simulator(BaseParamCollection())
    @df y plot!(plt, :t, :S, subplot = 1, color = :black, lw = 2)
    @df y plot!(plt, :t, :R, subplot = 2, color = :black, lw = 2)
    
    hyperZ = Truncated(Normal(1., 0.1), 0, Inf)
    y = replicates(DEBBase.simulator, BaseParamCollection(deb = DEBBaseParams(Z = hyperZ)), 100)
    @df y plot!(plt, :t, :S, group = :replicate, alpha = .25, subplot = 1, c = :viridis)
    @df y plot!(plt, :t, :R, group = :replicate, alpha = .25, subplot = 2, c = :viridis) 
    
    display(plt)
end

using BenchmarkTools
@benchmark [DEBBase.simulator(BaseParamCollection()) for _ in 1:100]

@benchmark replicates(DEBBase.simulator, BaseParamCollection(deb = DEBBaseParams(Z = hyperZ)), 1000)

macro replicates(expr::Expr, nreps::Int64)
    yhat = DataFrame()

    for replicate in 1:nreps
        yhat_i = eval(expr)
        yhat_i[!,:replicate] .= replicate
        append!(yhat, yhat_i)
    end

    return yhat
end








using DEBFigures
yhat = @replicates DEBBase.simulator(BaseParamCollection(deb = DEBBaseParams(Z = Truncated(Normal(1, 0.1), 0, Inf)))) 10

@df yhat plot(:t, :S, group = :replicate)

@df y plot(:t, :R, group = :replicate)

"""
DEBBase Agent. <br>
Each agent owns a reference to its associated parameter collection.
"""
mutable struct BaseAgent
    p::Base.RefValue{BaseParamCollection} # reference to the paramter collection
    u::ComponentVector
    du::ComponentVector

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

