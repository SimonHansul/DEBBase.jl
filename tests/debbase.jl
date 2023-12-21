using Revise 
using DEBBase
using Plots, StatsPlots

default(leg = false)


glb = GlobalBaseParams()
deb = DEBBaseParams()

out = simulator(glb, deb)

@less simulator(glb, deb) 

less(DEBBase.DEB!)

@df out plot(
    plot(:t, :S), 
    plot(:t, :R), 
    size = (600,350)
)


macro insert_function_body(func)
    body = func.body
    return quote
        $body
    end
end

@insert_function_body myfunction



@insert_function_body(DEBBase.run_model)

@enum drcfunct LL2 


deb.drc_funct_G
