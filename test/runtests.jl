using Pkg; Pkg.activate("test")

@info("Loading packjages")

begin
    using Revise 
    @time using DEBBase
    using SHUtils
    TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"")
    using Plots, StatsPlots, Plots.Measures
    default(titlefontsize = 10, lw = 1.5, leg = false)
    using Glob    
    using Plots, StatsPlots
    using DataFrames, DataFramesMeta
end

norm(x) = x ./ sum(x)

tests = glob("test/*.jl") |> 
x -> [splitpath(xi)[end] for xi in x] |>
x -> filter(f -> f != "runtests.jl", x)

for test in tests
    @info("Running $test")
    include(test)
end
