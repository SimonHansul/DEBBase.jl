using Pkg; Pkg.activate("test")

using Plots, StatsPlots
using Distributions
using DataFrames
default(leg = false)

using Revise
@time using DEBBase

defparams = BaseParamCollection()
defparams.deb.Z = Dirac(1.)
yhat = DEBBase.simulator(defparams)
@df yhat plot(:t, :S)

theta = BaseParamCollection()
theta.deb.Z = Truncated(Normal(1., 0.1), 0, Inf)
yhat = @replicates DEBBase.simulator(theta) 10

7
@df yhat plot(:t, :S, group = :replicate, color = 1)

