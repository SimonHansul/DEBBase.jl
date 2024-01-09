using Revise
using Plots, StatsPlots, Plots.Measures
default(titlefontsize = 10)

using DEBBase

plt = plot(
    layout = (1,3), leg = false, 
    title = ["Growth" "Reproduction" "Food density"], 
    leftmargin = 5mm, bottommargin = 5mm, 
    size = (800,350)
    )

for Xdot_in in [0.0125, 0.05, 0.1, 0.2, 0.125, .25, .5, 1.] .* 4800.
    out = simulator(
        BaseParamCollection(
            glb = GlobalBaseParams(Xdot_in = Xdot_in, t_max = 365.), 
            deb = DEBBaseParams(K_X = 12e3))
        )

    @df out plot!(plt, :t, :S, ylabel = "S", subplot = 1) 
    @df out plot!(plt, :t, :R, ylabel = "R", subplot = 2)
    @df out plot!(plt, :t, :X_p ./ GlobalBaseParams().V_patch, ylabel = "[X_p]", subplot = 3)
end
hline!(plt, [DEBBase.calc_S_max(DEBBaseParams())], linestyle = :dash, color = "gray", subplot = 1)
hline!(plt, [12e3], linestyle = :dash, color = "gray", subplot = 3)
display(plt)


Imax = 1.
S = -40.

I = S^(2/3) * Imax
S^(2/3)

(-40.0) ^(2/3)

