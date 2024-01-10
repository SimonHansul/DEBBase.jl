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
for Xdot_in in [0.0125, 0.05, 0.1, 0.2, 0.125, .25, .5, 1.] .* 4800.
    # generate the predidction
    out = simulator(
        BaseParamCollection(
            glb = GlobalBaseParams(Xdot_in = Xdot_in, t_max = 56.), 
            deb = DEBBaseParams(K_X = 12e3))
        )

    # plot the trajectories
    @df out plot!(plt, :t, :S, ylabel = "S", subplot = 1) 
    @df out plot!(plt, :t, :R, ylabel = "R", subplot = 2)
    @df out plot!(plt, :t, :X_p ./ GlobalBaseParams().V_patch, ylabel = "[X_p]", subplot = 3)
end
hline!(plt, [DEBBase.calc_S_max(DEBBaseParams())], linestyle = :dash, color = "gray", subplot = 1)
hline!(plt, [12e3], linestyle = :dash, color = "gray", subplot = 3)
display(plt)


using BenchmarkTools

@benchmark out = simulator(
        BaseParamCollection(
            glb = GlobalBaseParams(Xdot_in = 4800., t_max = 56.), 
            deb = DEBBaseParams(K_X = 12e3))
        )

#### nested sigmoid function

sig = DEBBase.sig
H_p = DEBBaseParams().H_p

function nsig(x, x_thr, y_steps) 
    return length(y_steps) - 1 + sig(x, x_thr[1], y_steps[1][1], y_steps[1][2]) + sum([-sig(x, x_thr[i+1], y_steps[i+1][2], y_steps[i+1][1]) for i in length(y_steps)-1])
end

plot(
    x -> 2 + sig(x, H_p, 2., 3.; β = 10.) - sig(x, 50., 2., 2.; β = 10.) - sig(x, 25., 1., 0.; β = 10.), 
    xlim = (0, 150.), xlabel = "H", ylabel = "life_stage", 
    title = "Nested sigmoid switch"
    )