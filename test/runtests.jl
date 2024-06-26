
if (abspath(PROGRAM_FILE) == @__FILE__) | occursin("terminalserver", abspath(PROGRAM_FILE)) 
    @info("Loading packages")
    using Pkg; Pkg.activate("test")
    using Plots, StatsPlots, Plots.Measures
    using StatsBase
    using SHUtils, Glob
    using DataFrames, DataFramesMeta
    using ProgressMeter
    default(leg = false, lw = 1.5)

    using Revise
    using DEBBase
end
TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"")

norm(x) = x ./ sum(x)
tests = glob("test/*.jl") |> 
x -> [splitpath(xi)[end] for xi in x] |>
x -> filter(f -> f != "runtests.jl", x)

include(tests[2])

y = simulator(DEBParamCollection())


for test in tests
    @info("Running $test")
    include(test)
end



using Revise
using DEBBase