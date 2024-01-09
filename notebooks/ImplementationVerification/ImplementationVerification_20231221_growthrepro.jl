using Revise
using Plots, StatsPlots
default(titlefontsize = 10)

using DEBBase

plt = plot(
    layout = (1,2), leg = false, 
    title = ["Growth" "Reproduction"], xlabel = "t", 
    size = (650,350)
    )

for Xdot_in in [.25, .5, 1.] .* 1200.
    out = simulator(
        BaseParamCollection(
            glb = GlobalBaseParams(Xdot_in = Xdot_in), 
            deb = DEBBaseParams(K_X = 12e3))
        )

    @df out plot!(plt, :t, :S, ylabel = "S", subplot = 1) 
    @df out plot!(plt, :t, :R, ylabel = "R", subplot = 2)
end

display(plt)
