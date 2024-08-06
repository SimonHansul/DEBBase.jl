begin
    @info("Loading packages")
    using Pkg; Pkg.activate("test")
    using Test
    using Plots, StatsPlots, Plots.Measures
    using DataFrames, DataFramesMeta
    using ProgressMeter
    default(leg = false, lw = 1.5)
    using OrdinaryDiffEq
    using Revise
    @time using DEBBase.DEBODE, DEBBase.Utils, DEBBase.Figures
    TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"")
end

using DEBBase
DEBBase.ABC.HierchPriors(
    hyperparams = [:Z],
    hyperpriors = [Truncated]
)


function f(x::Vector{Pair{Symbol,D}}) where D <: Distribution
    println(x)
end
