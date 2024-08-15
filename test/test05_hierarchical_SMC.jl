begin
    @info("Loading packages")
    using Pkg; Pkg.activate("test")
    using Test
    using Plots, StatsPlots, Plots.Measures
    using DataFrames, DataFramesMeta
    using ProgressMeter
    default(leg = false, lw = 1.5)
    using OrdinaryDiffEq
    using Distributions
    using Revise
    @time using DEBBase.ABC
    TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"")
end

priors = ABC.HierchPriors(
    :Z => Truncated(Normal(1, 1), 0, Inf),
    [:Z_1, :Z_2],
    [
        :Idot_max_rel => Truncated(Normal(1, 1), 0, Inf),
        :k_M => Truncated(Normal(0.1, 0.1), 0, Inf)
    ]
)



hierch_sample(priors)

using Distributions, Plots


hyperpriors = [Truncated(Normal(1, 1), 0, Inf)]



hypersample = rand.(hyperpriors)


plot(prior_cv)


Z = rand(Truncated(Normal(1, )))
