@info("Loading packages")
using Pkg; Pkg.activate("test")
using Plots, StatsPlots, Plots.Measures
<<<<<<< HEAD
default(leg = false, lw = 1.5)

using Test
using StatsBase
using DataFrames, DataFramesMeta, Distributions
using ProgressMeter
using BenchmarkTools
using Glob

using Revise
@time using DEBBase.DEBODE, DEBBase.ABC, DEBBase.Figures, DEBBase.Utils
=======
using StatsBase
using Glob
using DataFrames, DataFramesMeta
using ProgressMeter
default(leg = false, lw = 1.5)

using Revise
using DEBBase
using DEBBase.Figures

TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"")
>>>>>>> 041b237e17b3d5300ea06068f55701ba1bf5dfc2

norm(x) = x ./ sum(x)
tests = glob("test/*.jl") |> 
x -> [splitpath(xi)[end] for xi in x] |>
x -> filter(f -> f != "runtests.jl", x) |>
x -> filter(f -> occursin("test", f), x)

for test in tests
    @info("Running $test")
    include(test)
end

