using Revise 
@time using DEBBase
using Plots, StatsPlots, Plots.Measures

default(leg = false)

out = simulator(BaseParamCollection())

@df out plot(
    plot(:t, :S, ylabel = "S"), 
    plot(:t, :R, ylabel = "R"), 
    size = (600,350), 
    xlabel = "t", leftmargin = 5mm
)

