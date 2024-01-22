using Pkg; Pkg.activate("tests") # activate test environment
using Revise # load test packages
@time using DEBBase # should be added using "dev ." instead of "add .", so that current version is always loaded
using Plots, StatsPlots, Plots.Measures
default(titlefontsize = 10, lw = 1.5, leg = false)

using StatsPlots

p = BaseParamCollection()
out = simulator(p)

@df out plot(
    plot(:t, :S, ylabel = "S"), 
    plot(:t, :H, ylabel = "H"),
    plot(:t, :R, ylabel = "R"), 
    plot(:t, [:I, :I_p, :I_emb, :A], ylabel = "I"), 
    plot(:t, :life_stage, ylabel = "life_stage"),
    plot(:t, :X_p, ylabel = "X_p"), 
    plot(:t, :X_emb, ylabel = "X_emb"),
    xlabel = "t", 
    size = (700,450), 
    layout = (2,4)
)

norm(x) = x ./ (sum(x))

# prepare the plot
plt = plot(
    layout = grid(1,3, widths = norm([2, 1, 1])),
    leg = false, 
    title = ["Growth" "Rerpoduction" "Food density"], 
    leftmargin = 5mm, bottommargin = 6mm, 
    size = (1200,350), 
    xlabel = "Time (d)"
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
    @df out plot!(plt, :t, :S, ylabel = "S", subplot = 1, leg = :outertopleft, label = "Xdot_in = $(Xdot_in)") 
    @df out plot!(plt, :t, :R, ylabel = "R", subplot = 2)
    @df out plot!(plt, :t, :X_p ./ GlobalBaseParams().V_patch, ylabel = "[X_p]", subplot = 3, 
        yscale = :log10
        )
end
hline!(plt, [DEBBase.calc_S_max(DEBBaseParams())], linestyle = :dash, color = "gray", subplot = 1, label = "S_max")
display(plt)


using Dates
@time begin # forcing the code to run for 1 minute, to compare with output of @benchmark
    t0 = now()
    i = 0.
    while now() - t0 <= Second(60)
        out = simulator(
            BaseParamCollection(
                glb = GlobalBaseParams(Xdot_in = 4800., t_max = 21.), 
                deb = DEBBaseParams(K_X = 12e3))
            )
        i += 1
    end
    println("Carried out $i simulations in $Second((now() - t0)).")
end

using BenchmarkTools
@benchmark out = simulator(
        BaseParamCollection(
            glb = GlobalBaseParams(Xdot_in = 4800., t_max = 21.), 
            deb = DEBBaseParams(K_X = 12e3))
        )

using ProfileView
ProfileView.@profview(simulator(BaseParamCollection()))
