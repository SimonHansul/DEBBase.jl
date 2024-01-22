using Pkg; Pkg.activate("tests")
using Revise
using DEBBase
using Plots, StatsPlots

p = BaseParamCollection()
isolate_pmoas!(p.deb, ["A"])
p.glb.C_W = [0.5]
p.deb.k_D_A = [1.]
p.deb.drc_params_A = [(1., 2.)]
simout = simulator(p)

@df simout plot(
    plot(:t, :S),
    plot(:t, :R),
    plot(:t, :D_A_1)
)
