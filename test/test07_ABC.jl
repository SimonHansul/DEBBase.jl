using Pkg; Pkg.activate("test")

using Plots, StatsPlots, StatsBase
default(legendtitlefontsize = 10)
using DataFrames, DataFramesMeta, CSV, YAML

using Revise

using DEBBase.Utils
using DEBBase.DEBODE
using DEBBase.ABC
using DEBBase.Figures


#=
## Loading a data config file
=#

data = Utils.data_from_config("test/config/data_config_example.yml") 

function plot_data(data::Utils.Dataset; kwargs...)

    plt = plot(layout = (1, 2); kwargs...)
    
    @df data.time_resolved["growth_agg"] plot!(plt, subplot = 1,
        :t_birth, :drymass_mean, group = :food_level, 
        yerr = :drymass_sd,
        c = :black, lw = 2, marker = [:diamond :circle], linestyle = [:dash :solid],
        xlabel = "Time since birth (d)", ylabel = "Mean dry mass (mg)", legendtitle = "Food level (mg/d)"
        )

    @df data.time_resolved["repro_agg"] plot!(plt, subplot = 2,
        :t_birth, :cum_repro_mean, group = :food_level, 
        yerr = :cum_repro_sd,
        c = :black, lw = 2, marker = [:diamond :circle], linestyle = [:dash :solid],
        xlabel = "Time since birth (d)", ylabel = "Mean cumulative \n reproduction (#)", 
        legendtitle = "Food level (mg/d)"
        )
    
    return plt
end

plot_data(data, size = (800,500))

#=
## Defining an initial guess
=#

begin
    Z = maximum(data.time_resolved["growth_agg"].drymass_mean) / DEBODE.calc_S_max(SpeciesParams())

    intguess = Params()
    spc = intguess.spc
    spc.Idot_max_rel_0 *= Z^(1/3)#
    spc.Idot_max_rel_emb_0 = spc.Idot_max_rel_0
    spc.X_emb_int_0 *= Z
    spc.H_p_0 = 0.03

    yhat = DEBODE.simulator(intguess)

    plot_data(data)

    @df yhat plot!(:t, :S, lw = 2, label = "Initial guess")
    @df yhat plot!(:t, :R ./ spc.X_emb_int_0, lw = 2, subplot = 2, label = "")
end




#=
## Defining the simulator
=#

using DEBBase.ParamStructs

function simulate_data(
    params::AbstractParamCollection, 
    parnames, 
    parvals; 
    food_levels = [0.1, 0.4]
    )

    # assigning parameter samples
    for (par,val) in zip(parnames,parvals)
        eval(:(params.spc.$par = $val))
    end

    sim = DataFrame()
    for food_level in food_levels 
        params.glb.Xdot_in = food_level
        sim_i = DEBODE.simulator(params)

        sim_i[!,:food_level] .= food_level 
        sim_i[!,:t_birth] = sim_i.t .- DEBODE.age_at_birth(sim_i) 

        append!(sim, sim_i)
    end
    

    yhat = Dataset()
    yhat.time_resolved["growth_agg"]

    return yhat
end


data.time_resolved["growth_agg"]


yhat = simulate_data(intguess, [], [])

@df yhat plot(:t, :S, group = :food_level)





