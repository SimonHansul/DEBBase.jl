

using Pkg
Pkg.activate("./test")
Pkg.add("")

using Revise
using StatsPlots
@time using DEBBase
using Plots, StatsPlots
default(leg = false)


# lower β in combination with sig?

p = BaseParamCollection()
p.glb.t_max = 21.
out = simulator(p)

@df out plot(
    plot(:t, :S),
    plot(:t, :X_emb),
    plot(:t, :X_p),
    plot(:t, [:I, :I_emb, :I_p])
)


using Cthulhu
using BenchmarkTools

using DEBBase
simulator(p)

@benchmark simulator(p)

