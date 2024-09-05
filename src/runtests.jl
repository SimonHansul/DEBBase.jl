
@info("Loading packages")
using Pkg; Pkg.activate("test")
using Plots, StatsPlots, Plots.Measures
default(leg = false, lw = 1.5)

using StatsBase
using DataFrames, DataFramesMeta
using ProgressMeter

using Glob

using Revise
using DEBBase.ODE, DEBBase.ABC, DEBBase.Figures, DEBBase.Utils

TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"")
println(TAG)

norm(x) = x ./ sum(x)
tests = glob("test/*.jl") |> 
x -> [splitpath(xi)[end] for xi in x] |>
x -> filter(f -> f != "runtests.jl", x)


for test in tests
    @info("Running $test")
    include(test)
end

