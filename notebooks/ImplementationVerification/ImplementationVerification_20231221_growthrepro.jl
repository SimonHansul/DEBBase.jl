using Revise
using Plots, StatsPlots, Plots.Measures
default(titlefontsize = 10)

default(leg = false)

@time using DEBBase

# prepare the plot
plt = plot(
    layout = (1,3), leg = false, 
    title = ["Growth" "Reproduction" "Food density"], 
    leftmargin = 5mm, bottommargin = 5mm, 
    size = (800,350)
    )

# iterate over nutrient input concentrations
Xdot_in = 4800.
for _ in 1:5
    Xdot_in /= 2
    # generate the predidction
    out = simulator(
        BaseParamCollection(
            glb = GlobalBaseParams(Xdot_in = Xdot_in, t_max = 56.), 
            deb = DEBBaseParams(K_X = 12e3))
        )

    # plot the trajectories
    @df out plot!(plt, :t, :S, ylabel = "S", subplot = 1, leg = true, label = Xdot_in) 
    @df out plot!(plt, :t, :R, ylabel = "R", subplot = 2)
    @df out plot!(
        plt, :t, :X_p ./ GlobalBaseParams().V_patch, ylabel = "[X_p]", subplot = 3, 
        yscale = :log10
        )
end
hline!(plt, [DEBBase.calc_S_max(DEBBaseParams())], linestyle = :dash, color = "gray", subplot = 1, label = "S_max")
hline!(plt, [12e3], linestyle = :dash, color = "gray", subplot = 3)
display(plt)


using BenchmarkTools
using DEBBase

@benchmark out = simulator(
        BaseParamCollection(
            glb = GlobalBaseParams(Xdot_in = 4800., t_max = 56.), 
            deb = DEBBaseParams(K_X = 12e3))
        )

@time out = simulator(
    BaseParamCollection(
        glb = GlobalBaseParams(Xdot_in = 4800., t_max = 56.), 
        deb = DEBBaseParams(K_X = 12e3))
    )

@df out plot(
    plot(:t, :S), 
    plot(:t, :R)
)




plot(
    x -> 2 + sig(x, H_p, 2., 3.; β = 10.) - sig(x, 50., 2., 2.; β = 10.) - sig(x, 25., 1., 0.; β = 10.), 
    xlim = (0, 150.), xlabel = "H", ylabel = "life_stage", 
    title = "Nested sigmoid switch"
    )