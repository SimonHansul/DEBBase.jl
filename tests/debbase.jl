using Revise 
@time using DEBBase
using Plots, StatsPlots, Plots.Measures

default(leg = false)
glb = GlobalBaseParams()
deb = DEBBaseParams()
out = simulator(glb, deb)

@df out plot(
    plot(:t, :S, ylabel = "S"), 
    plot(:t, :R, ylabel = "R"), 
    size = (600,350), 
    xlabel = "t", leftmargin = 5mm
)

a, b, c = [1, 2, 3, 4]

a, b, c