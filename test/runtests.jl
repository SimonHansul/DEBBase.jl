@info("Loading packages")
using Pkg; Pkg.activate("test")
using Plots, StatsPlots, Plots.Measures
default(leg = false, lw = 1.5)

using Test
using StatsBase
using DataFrames, DataFramesMeta, Distributions
using ProgressMeter
using BenchmarkTools
using Glob

using Revise
@time using DEBBase.DEBODE, DEBBase.ABC, DEBBase.Figures, DEBBase.Utils

yhat = simulator(DEBParamCollection(), saveat = 1/24)
@df yhat plot(
    plot(:t, [:embryo :juvenile :adult]),
    plot(:t, :X_emb), 
    plot(:t, :S), 
    plot(:t, :H),
    plot(:t, :A)
)


norm(x) = x ./ sum(x)
tests = glob("test/*.jl") |> 
x -> [splitpath(xi)[end] for xi in x] |>
x -> filter(f -> f != "runtests.jl", x) |>
x -> filter(f -> occursin("test", f), x)

for test in tests
    @info("Running $test")
    include(test)
end

