using Pkg; Pkg.activate("tests")
using Plots, Plots.Measures, StatsPlots
using DataFrames, DataFramesMeta
using BenchmarkTools

using Revise
@time using DEBBase

default(leg = false, lw = 1.5)

p = BaseParamCollection()
simout = simulator(p)

plt = @df simout plot(
    plot(:t, :S, ylabel = "S [μg C]"), 
    plot(:t, :R ./ p.deb.X_emb_int, ylabel = "R [#]"), 
    xlabel = "t", color = :black
)
savefig(plt, "plots/basetest_growthrepo.png")
@benchmark simulator(p)