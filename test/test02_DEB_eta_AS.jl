using Pkg; Pkg.activate("tests")
using Plots, StatsPlots, Plots.Measures
default(titlefontsize = 10, lw = 1.5, leg = false)
using Revise 
@time using DEBBase
using DataMonk
const TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"")

norm(x) = x ./ (sum(x))

begin # effect of food input
    plt = plot( # prepare the plot
        layout = grid(1,3, widths = norm([2, 1, 1])),
        leg = false, 
        title = ["Growth" "Reproduction" "Food density"], 
        leftmargin = 5mm, bottommargin = 6mm, 
        size = (1200,350), 
        xlabel = "Time (d)"
        )

    # iterate over nutrient input concentrations
    let eta_AS_vec = collect(0.1:0.1:1)
        for eta_AS in eta_AS_vec
            # generate the predidction
            out = simulator(
                BaseParamCollection(
                    glb = GlobalBaseParams(t_max = 56.), 
                    deb = DEBBaseParams(eta_AS = eta_AS))
                )

            # plot the trajectories
            @df out plot!(plt, :t, :S, ylabel = "S", subplot = 1, leg = :outertopleft, label = "eta_AS = $(eta_AS)") 
            @df out plot!(plt, :t, :R, ylabel = "R", subplot = 2)
            @df out plot!(plt, :t, :X_p ./ GlobalBaseParams().V_patch, ylabel = "[X_p]", subplot = 3, 
                yscale = :log10
                )
        end
        hline!(plt, [DEBBase.calc_S_max(DEBBaseParams())], linestyle = :dash, color = "gray", subplot = 1, label = "S_max")
        display(plt)
        savefig(plt, "plots/$(TAG).png")
    end
end

