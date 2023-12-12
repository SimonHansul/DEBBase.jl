using Revise
@time using DEBBase
using Plots, Plots.Measures, StatsPlots
using DataFrames
using BenchmarkTools

default(leg = false, lw = 1.5)

glb = GlobalBaseParams()
anm = DEBBaseParams()
simout = DEBBase.run_model(glb, anm)

