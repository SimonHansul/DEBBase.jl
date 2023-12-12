using Infiltrator
using Revise
@time using DEBBase
using Plots, Plots.Measures

default(leg = false, lw = 1.5)

glb = GlobalBaseParams()
anm = DEBBaseParams()
sol = DEBBase.run_model(glb, anm)

out = hcat(sol.u...)'
plot!(
    plot(sol.t, out[:,1], ylabel = "X"),
    plot(sol.t, out[:,2], ylabel = "X_emb"),
    plot(sol.t, out[:,3], ylabel = "S"),
    plot(sol.t, out[:,4], ylabel = "H"),
    plot(sol.t, out[:,5] ./ anm.X_emb_int, ylabel = "R"), 
    size = (800,450), leftmargin = 5mm
)
hline!([DEBBase.calc_S_max(anm)], subplot = 3, linestyle = :dash, color = "gray")



