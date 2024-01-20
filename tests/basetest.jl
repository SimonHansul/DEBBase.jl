using Pkg; Pkg.activate("tests")
using Plots, Plots.Measures, StatsPlots
using BenchmarkTools
default(leg = false, lw = 1.5)

using Revise
@time using DEBBase

p = BaseParamCollection()
p.deb.k_D_G = [0.]
p.deb.k_D_M = [0.]
p.deb.k_D_A = [1.]
p.deb.k_D_R = [0.]
p.deb.drc_params_A = [(1., 2.)]

p.glb.C_W = [0.]
simout = simulator(p)
plt = @df simout plot(
    plot(:t, :S, ylabel = "S [μg C]"), 
    plot(:t, :R ./ p.deb.X_emb_int, ylabel = "R [#]"), 
    xlabel = "t", color = :black
)

p.glb.C_W = [0.5]
simout = simulator(p)
@df simout plot!(:t, :S, subplot = 1)

@df simout plot(:t, :D_A_1, ylim = (0, Inf))