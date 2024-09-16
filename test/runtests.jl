@info("Loading packages")
using Pkg; Pkg.activate("test")
using Plots, StatsPlots, Plots.Measures
using StatsBase
using Glob
using DataFrames, DataFramesMeta
using ProgressMeter
default(leg = false, lw = 1.5)

using Revise
using DEBBase
using DEBBase.Figures

TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"")

norm(x) = x ./ sum(x)
tests = glob("test/*.jl") |> 
x -> [splitpath(xi)[end] for xi in x] |>
x -> filter(f -> f != "runtests.jl", x)

for test in tests
    @info("Running $test")
    include(test)
end



