using Revise 
using DEBBase
using Plots, StatsPlots

default(leg = false)


glb = GlobalBaseParams()
deb = DEBBaseParams()

out = DEBBase.run_model(glb, deb)

@df out plot(
    plot(:t, :S), 
    plot(:t, :R), 
    size = (600,350)
)

@enum drcfunct LL2 


deb.drc_funct_G
