using Pkg; Pkg.activate("tests")
using Plots, StatsPlots
using BenchmarkTools
using Revise
begin
    @time using DEBBase

    p = BaseParamCollection(glb = GlobalBaseParams(C_W = [1.]));
    y = simulator(p);
    @benchmark simulator(p)
end
isolate_pmoas!(p.deb, ["G"])
p.glb.C_W = [0.5]
p.deb.k_D_A = [1.]
p.deb.drc_params_A = [(1., 2.)]
simout = simulator(p)

@df simout plot(
    plot(:t, :S),
    plot(:t, :R),
    plot(:t, :D_A_1)
)
