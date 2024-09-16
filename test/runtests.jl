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

using Revise
@time using DEBBase.DEBODE, DEBBase.Utils, DEBBase.Figures

norm(x) = x ./ sum(x)
tests = glob("test/*.jl") |> 
x -> [splitpath(xi)[end] for xi in x] |>
<<<<<<< HEAD
x -> filter(f -> f != "runtests.jl", x)
=======
x -> filter(f -> f != "runtests.jl", x) |>
x -> filter(f -> occursin("test", f), x)
>>>>>>> 0f934f80debb4f2780b085f08fbe9d06af5d466e

for test in tests
    @info("Running $test")
    include(test)
end

SpeciesParams()


