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

# For easier access while developing this code, we can as well create aliases for some of the entries 
# in the dataset. 
# This does by the way not allocate new memory.

growth = data.time_resolved["growth_agg"] 
repro = data.time_resolved["repro_agg"]
stats_growth = data.scalar["growth_stats_agg"]
stats_repro = data.scalar["repro_stats_agg"]
traits = data.scalar["misc_traits"]

# We also define a function to plot the time-resolved data.

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

We define an initial guess, starting with the default parameters. 
Based on the observed maximum drymass, we calculate a zoom factor and apply it to the default parameters to obtain our initial guess. 

Comparing the prediction of the initial guess with the observed values, we should be in the right order of magnitude.
=#

begin
    Z = maximum(data.time_resolved["growth_agg"].drymass_mean) / DEBODE.calc_S_max(SpeciesParams())

    intguess = Utils.params_from_config(Params, "test/config/param_config_abctest.yml") # loads a config file containing the global settings. we could as well enter the values here, but the config file can help with reproducibility

    spc = intguess.spc
    spc.Idot_max_rel_0 *= Z^(1/3)#
    spc.Idot_max_rel_emb_0 = spc.Idot_max_rel_0
    spc.X_emb_int_0 *= Z
    spc.H_p_0 = 0.03
    spc.K_X_0 = 1.0 # guess for the half-saturation constant (mg/L)

    yhat = DEBODE.simulator(intguess; reltol = 1e-10)

    plot_data(data)

    @df yhat plot!(:t, :S, lw = 2, label = "Initial guess")
    @df yhat plot!(:t, :R ./ spc.X_emb_int_0, lw = 2, subplot = 2, label = "")
end

#=
## Defining the simulator

The `simulate_data` function first generates the raw simulation output for both food levels based on a parameter sample. <br>
But this is not enough, because our data contains values which are not directly comparable with the raw simulation output. <br<
The second portion of `simulate_data` therefore processes the simulation output, and finally returns it as a `Dataset` object.
=#

using DEBBase.ParamStructs
using DataStructures

 function simulate_data(
    defaultparams::AbstractParamCollection, 
    parnames::Vector{Symbol}, # names of estimated parameters
    parvals::Vector{R}; # samples parameter values
    food_levels = [0.1, 0.4], # simulated food levels
    return_raw = false, # optionally return the raw simulation output
    ) where R <: Real

    params = deepcopy(defaultparams)

    # assigning parameter samples
    for (par,val) in zip(parnames,parvals)
        setproperty!(params.spc, par, val)
        #eval(:(params.spc.$par = $val))
    end

    params.spc.k_J_0 = ((1 - params.spc.kappa_0) / params.spc.kappa_0) * params.spc.k_M_0
    # generating the simulation output
    sim = DataFrame()
    for food_level in food_levels 
        params.glb.Xdot_in = food_level
        sim_i = DEBODE.simulator(params, saveat = 1)

        sim_i[!,:food_level] .= food_level 
        sim_i[!,:t_birth] = sim_i.t .- round(DEBODE.age_at_birth(sim_i)) 

        append!(sim, sim_i)
    end

    if return_raw
        return sim
    end

    # processing the simulation output

    yhat = Utils.Dataset() # an empty dataset

    # adding simulated time-resolved data

    yhat.time_resolved["growth_agg"] = @chain sim begin
        @select(:t_birth, :food_level, :S, :R)
        @transform(:drymass_mean = @. :S)
        unique
    end

    yhat.time_resolved["repro_agg"] = @chain sim begin
        @select(:t, :t_birth, :food_level, :R)
        @transform(:cum_repro_mean = @. trunc(:R / params.spc.X_emb_int_0)) # NOTE: this is only valid if we don't simulate individual variability on X_emb
        @transform(:t_birth = :t) # instead of shifting the repro values, we shift the time-values - since t_birth = t - age_at_birth and we shift by age_at_birth, we can just say t_birth = t
        select(Not(:t)) # not needed anymore
    end

    # adding simulated scalar data

    yhat.scalar["growth_stats_agg"] = @chain yhat.time_resolved["growth_agg"] begin
        groupby(:food_level)
        combine(_) do df
            return DataFrame(
                max_drymass_mean = maximum(df.drymass_mean),
                drymass_at_birth_mean =  @subset(df, ismin(df.t_birth)).drymass_mean[1]
            )
        end
    end

    yhat.scalar["repro_stats_agg"] = @chain yhat.time_resolved["repro_agg"] begin
        groupby(:food_level)
        combine(_) do df
            return DataFrame(
                birth_to_first_hatch_mean = @subset(df, :cum_repro_mean .> 0).t_birth |> Utils.robustmin,
                final_cum_repro_mean = @subset(df, :t_birth .== maximum(sim.t_birth)).cum_repro_mean |> x-> length(x)>0 ? x[1] : NaN
            )
        end
    end

    yhat.scalar["misc_traits"] = OrderedDict(
        "age_at_birth" => Dict("value" => DEBODE.age_at_birth(@subset(sim, :food_level .== maximum(:food_level)))),
        "egg_drymass" => Dict("value" => params.spc.X_emb_int_0)
    )

    return yhat
end

simulate_data(params::AbstractParamCollection; kwargs...) = simulate_data(params, Symbol[], Real[]; kwargs...)

using Test
@testset begin
    @info "Simulating a Dataset"

    global sim = simulate_data(intguess; return_raw = true) # this would return the raw ODE output
    global yhat = simulate_data(intguess) # this returns the simulated dataset

    @test sim isa DataFrame
    @test yhat isa ABC.AbstractDataset
end


#=
## Defining the distance function

The symmetric bounded loss function is used as 
    our "inner" distance function. 
That means, the symmetric bounded loss is applied to each part of the dataset, 
then combined to calculate the total distance.
=#

using Test

observed = data
predicted = yhat

@testset begin
    @info "Computing loss for time-resolved data"

    loss = ABC.compute_time_resolved_loss("repro_agg", predicted, observed)

    @test isfinite(loss)
end

predicted.time_resolved["repro_agg"]

@testset begin
    @info "Computing loss for scalar data in tabular format"
    name = "growth_stats_agg"

    observed_df = observed.scalar[name]
    predicted_df = predicted.scalar[name]
    grouping_vars = observed.grouping_vars["scalar"][name]
    response_vars = observed.response_vars["scalar"][name]

    loss = ABC.compute_scalar_loss(predicted_df, observed_df, response_vars, grouping_vars)

    @test isfinite(loss)
end

@testset begin
    @info "Computing loss for scalar data in dict format"
    name = "misc_traits"

    observed_dict = observed.scalar[name]
    predicted_dict = predicted.scalar[name]
    grouping_vars = observed.grouping_vars["scalar"][name]
    response_vars = observed.response_vars["scalar"][name]

    println(grouping_vars)

    loss = ABC.compute_scalar_loss(
        predicted_dict, observed_dict,
        response_vars, grouping_vars
        )

    @test isfinite(loss)
end

@testset begin
    @info "Computing total loss"

    loss = ABC.compute_loss(predicted, observed)
    @test isfinite(loss)
end


#=
Below we define the priors as truncated Normal distributions with a constant CV of 1,000% around the initial guess.
=#

begin 

    fitted_params = [
        #:X_emb_int_0,
        #:Idot_max_rel_emb_0,
        :K_X_0,
        :Idot_max_rel_0,
        :kappa_0,
        :k_M_0,
        :eta_AS_0,
        :H_p_0
    ]

    cvs = fill(1.0, length(fitted_params))
    #cvs[1] = 0.1
    #cvs[2] = 0.5

    upper_limits = [
        #Inf,
        #Inf, 
        Inf, 
        Inf, 
        1,
        Inf,
        1,
        Inf
    ]

    priors = Priors(
        fitted_params, 
        [
            deftruncnorm(eval(:(intguess.spc.$param)), cv, u = upper) for (param,upper,cv) in zip(fitted_params,upper_limits,cvs)
        ]
    )


    plot(priors, layout = (3, 3), size = (800,500), leg = false)
end

#=
Before we scale up the parameter inference, 
we can try a small SMC run with a strict rejection threshold and fewer samples. 

This will give us an indiciation wether the loss function 
=#

begin

@time smc = SMC(
    priors,
    intguess,
    simulate_data,
    ABC.compute_loss,
    data;
    n_pop = 5_000, 
    q_eps = 0.1,
    k_max = 0
)

bestfit = deepcopy(intguess)
ABC.posterior_sample!(
    bestfit.spc, 
    @subset(smc.accepted, :distance .== minimum(:distance))
)

yhat_bestfit = simulate_data(bestfit)
sim_bestfit = simulate_data(bestfit, return_raw = true)
plt_bestfit = plot_data(data)
@df yhat_bestfit.time_resolved["growth_agg"] plot!(subplot = 1, :t_birth, :drymass_mean, group = :food_level)
@df yhat_bestfit.time_resolved["repro_agg"] plot!(subplot = 2, :t_birth, :cum_repro_mean, group = :food_level)
display(plt_bestfit)

end

data.scalar["growth_stats_agg"]
yhat_bestfit.scalar["growth_stats_agg"]

@df sim plot(:t_birth, :f_X, group = :food_level)