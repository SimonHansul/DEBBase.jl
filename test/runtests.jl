@time begin 
    @info("Loading packages")
    using Pkg; Pkg.activate("test")
    using Plots, StatsPlots, Plots.Measures
    using StatsBase
    using Glob
    using DataFrames, DataFramesMeta
    using ProgressMeter
    using Distributions
    default(leg = false, lw = 1.5)
    using Test
end
using Revise
using DEBBase.DEBODE, DEBBase.Utils
using DEBBase.Figures


TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"")

norm(x) = x ./ sum(x)
tests = glob("test/*.jl") |> 
x -> [splitpath(xi)[end] for xi in x] |>
x -> filter(f -> f != "runtests.jl", x) |>
x -> filter(f -> occursin("test", f), x)

# show debugging statements
# set to "" to disable
ENV["JULIA_DEBUG"] = Main 

for test in tests
    @info("Running $test")
    include(test)
end

