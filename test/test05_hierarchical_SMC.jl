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
