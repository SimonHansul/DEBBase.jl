#=
## Estimating DEB parameters from control data at a single food level
=#

#=
## Setup
=#

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

# We only use the highest food level for this test to simplify the test a little

data.time_resolved["growth_agg"] = @subset(data.time_resolved["growth_agg"], :food_level .== 0.4)
data.time_resolved["repro_agg"] = @subset(data.time_resolved["repro_agg"], :food_level .== 0.4)

data.scalar["growth_stats_agg"] = @subset(data.scalar["growth_stats_agg"], :food_level .== 0.4)
data.scalar["repro_stats_agg"] = @subset(data.scalar["repro_stats_agg"], :food_level .== 0.4)

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

# TODO: this can be done more generically, so we don't have to re-write the entire simulator every time
function compute_growth_agg(sim::DataFrame)
    return @chain sim begin
        @select(:t_birth, :food_level, :S, :R)
        @transform(:drymass_mean = @. :S)
        unique
    end
end

function compute_repro_agg(sim::DataFrame, params::AbstractParamCollection)
    return @chain sim begin
        @select(:t, :t_birth, :food_level, :R)
        @transform(:cum_repro_mean = @. trunc(:R / params.spc.X_emb_int_0)) # NOTE: this is only valid if we don't simulate individual variability on X_emb
        @transform(:t_birth = :t) # instead of shifting the repro values, we shift the time-values - since t_birth = t - age_at_birth and we shift by age_at_birth, we can just say t_birth = t
        select(Not(:t)) # not needed anymore
    end
end

function compute_growth_stats_agg(yhat::Utils.AbstractDataset)
    @chain yhat.time_resolved["growth_agg"] begin
        groupby(:food_level)
        combine(_) do df
            return DataFrame(
                max_drymass_mean = maximum(df.drymass_mean),
                drymass_at_birth_mean =  @subset(df, ismin(df.t_birth)).drymass_mean[1]
            )
        end
    end
end

function compute_repro_stats_agg(yhat::Utils.AbstractDataset)
    return @chain yhat.time_resolved["repro_agg"] begin
        groupby(:food_level)
        combine(_) do df
            return DataFrame(
                birth_to_first_hatch_mean = @subset(df, :cum_repro_mean .> 0).t_birth |> Utils.robustmin,
                final_cum_repro_mean = @subset(df, :t_birth .== maximum(sim.t_birth)).cum_repro_mean |> x-> length(x)>0 ? x[1] : NaN
            )
        end
    end
end

function compute_misc_traits(sim::DataFrame, params::AbstractParamCollection)
    let S_max_theo = DEBODE.calc_S_max(params.spc)
        return OrderedDict(
            "age_at_birth" => Dict("value" => DEBODE.age_at_birth(@subset(sim, :food_level .== maximum(:food_level)))),
            "egg_drymass" => Dict("value" => maximum(sim.X_emb)),
            "S_max_rel" => Dict("value" => maximum(sim.S)/S_max_theo)
        )
    end
end

 function simulate_data(
    defaultparams::AbstractParamCollection, 
    parnames::Vector{Symbol}, # names of estimated parameters
    parvals::Vector{R}; # samples parameter values
    food_levels = [0.1, 0.4], # simulated food levels
    return_raw = false, # optionally return the raw simulation output
    ) where R <: Real

    @assert !(:Idot_max_rel_emb_0 in parnames) "Idot_max_rel_emb_0 is provided as parameter, but this value is calculated internally"

    params = deepcopy(defaultparams) # this is unfortunately still necessary if we use multithreading - might be easier to replace if we eventually switch to distributed computing

    # assigning parameter samples
    for (par,val) in zip(parnames,parvals)
        setproperty!(params.spc, par, val)
        #eval(:(params.spc.$par = $val))
    end

    params.spc.k_J_0 = ((1 - params.spc.kappa_0) / params.spc.kappa_0) * params.spc.k_M_0
    params.spc.Idot_max_rel_emb_0 = params.spc.Idot_max_rel_0

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

    yhat.time_resolved["growth_agg"] = compute_growth_agg(sim)
    yhat.time_resolved["repro_agg"] = compute_repro_agg(sim, params)

    # adding simulated scalar data

    yhat.scalar["growth_stats_agg"] = compute_growth_stats_agg(yhat)
    yhat.scalar["repro_stats_agg"] = compute_repro_stats_agg(yhat)

    yhat.scalar["misc_traits"] = compute_misc_traits(sim, params)

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

observed = data
predicted = yhat

@testset begin
    @info "Computing loss for time-resolved data"

    loss = ABC.compute_time_resolved_loss("repro_agg", predicted, observed)

    @test isfinite(loss)
end

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
        :X_emb_int_0,
        :Idot_max_rel_0,
        #:Idot_max_rel_emb_0,
        #:K_X_0,
        :kappa_0,
        :eta_AS_0,
        :k_M_0,
        :H_p_0
    ]

    cvs = fill(1.0, length(fitted_params))

    cvs[1] = 1.0 # we have a pretty good idea of the egg from data - narrowing this prior
    cvs[2] = 1.0 # same for maximum ingestion

    lower_limits = fill(0., length(fitted_params))
    upper_limits = [
        Inf,
        Inf, 
        #Inf, 
        #Inf,
        1,
        1,
        Inf,
        Inf
    ]

    priors = Priors(
        fitted_params, 
        [
            deflognorm(
                eval(:(intguess.spc.$param)), 
                sigma, 
                u = upper
                ) 
                for (param,upper,sigma) in zip(fitted_params,upper_limits,cvs)
        ]
    )


    plot(priors, layout = (3, 3), size = (800,500), leg = false, ylabel = gridylabel("Density", 3, 3))
end


#=
## Prior predictive check
=#


using DEBBase.ABC

prior_predictions = ABC.prior_predictice_check(
    priors, 
    intguess, 
    simulate_data, 
    data, 
    ABC.compute_loss
    );

begin
    prior_predictive_plot = plot_data(data)

    for p in prior_predictions.predictions
        @df p.time_resolved["growth_agg"] plot!(
            prior_predictive_plot, subplot = 1,
            :t_birth, :drymass_mean, group = :food_level, 
            color = :gray, label = ""
        )

        @df p.time_resolved["repro_agg"] plot!(
            prior_predictive_plot, subplot = 2,
            :t_birth, :cum_repro_mean, group = :food_level,
            color = :gray, label = ""
        )

    end

    prior_predictive_plot
end

plot(
    histogram(prior_predictions.distances, leg = false, fill = :gray, fillalpha = .5, lw = 1.5), 
    xlabel = "Loss", ylabel = "Count", 
    title = "Distribution of losses in prior predictive check \n 
        $(round(100*sum(isfinite.(prior_predictions.distances)/length(prior_predictions.distances))))% valid losses",
    titlefontsize = 10
)

#=
## Parameter inference

Before we scale up the parameter inference, 
we can try a small SMC run with a strict rejection threshold and fewer samples. 

=#

begin
    @time "Inferring posteriors using SMC" smc = SMC(
        priors,
        intguess,
        simulate_data,
        ABC.compute_loss,
        data;
        n_pop = 7_500, 
        q_eps = .2,
        k_max = 5
    )

    # quick check of the point estimate

    bestfit = deepcopy(intguess)
    ABC.posterior_sample!(
        bestfit.spc, 
        @subset(smc.accepted, :distance .== minimum(:distance))
    )

    yhat_bestfit = simulate_data(bestfit)
    sim_bestfit = simulate_data(bestfit; return_raw = true)

    plt_bestfit = plot_data(data, size = (800,500), leg = [false :topleft])
    @df yhat_bestfit.time_resolved["growth_agg"] plot!(subplot = 1, :t_birth, :drymass_mean, group = :food_level)
    @df yhat_bestfit.time_resolved["repro_agg"] plot!(subplot = 2, :t_birth, :cum_repro_mean, group = :food_level)
    display(plt_bestfit)
end

using Plots.Measures
begin
    distplot = plot(
        smc.priors, 
        size = (800,400), 
        label = "prior", 
        bottommargin= 5mm, 
        layout = (2,3), 
        ylabel = gridylabel("Density", 2, 3), 
        leftmargin = 5mm
        )
    for (i,param) in enumerate(priors.params)
        density!(
            smc.accepted[:,param], 
            weights = Weights(smc.accepted.weight),
            color = :gray, lw = 1.5, fill = true, fillalpha = .25,
            label = "posterior",
            subplot = i
            )
    end
    distplot        
end

post_pred = ABC.posterior_predictions(
    intguess,
    simulate_data,
    smc.accepted,
    priors
)

begin    
    vpcplot = plot_data(data)
    
    for yhat in post_pred
        @df yhat.time_resolved["growth_agg"] plot!(subplot = 1, :t_birth, :drymass_mean, group = :food_level, alpha = .05, color = 1, label = "")
        @df yhat.time_resolved["repro_agg"] plot!(subplot = 2, :t_birth, :cum_repro_mean, group = :food_level, alpha = .05, color = 1, label = "")
    end

    vpcplot
end


#=
## Plausability check
=#


sim_bestfit = combine(groupby(sim_bestfit, :food_level)) do df
    df[!,:dI] = diffvec(df.I)
    return df
end


@df sim_bestfit plot(
    plot(:t, :X_p, group = :food_level),
    plot(:t, :f_X, group = :food_level, ylim = (0, 1.01)),
    plot(:t, :dI, group = :food_level), 
    scatter(:X_p ./ intguess.glb.V_patch, :f_X, xlim = (-1, 10)), 
    lw = 2, 
    ylabel = ["Xₚ" "fₓ" "dI" "fₓ"], 
    xlabel = ["t" "t" "t" "[Xₚ]"]
)


@df yhat_bestfit.time_resolved["growth_agg"] plot(
    :t_birth, :drymass_mean, group = :food_level
    )



#=

## Model fitting with Nelder Mead

Below is an example that uses the Optim package and Nelder-Mead.
=#

using Optim

begin # adjusted weights to only fit to growth + misc traits
    f(x) = simulate_data(intguess, priors.params, x) |>     # minimization function
    x -> ABC.compute_loss(x, data)

    x0 = [mode(p.untruncated) for p in priors.priors]    # initial guessses
    @time "Fitting model using Nelder Mead" optim = optimize(  # performing the optimization
        f, x0, NelderMead(maxiter = 100)
        )   
    bestfit = optim.minimizer   # retrieving the estimates
    yhat_bestfit = simulate_data(intguess, priors.params, bestfit)   # plotting the prediction

    plt_bestfit = plot_data(data, size = (800, 350))
    @df yhat_bestfit.time_resolved["growth_agg"] plot!(subplot = 1, :t_birth, :drymass_mean, group = :food_level, lw = 2)
    @df yhat_bestfit.time_resolved["repro_agg"] plot!(subplot = 2, :t_birth, :cum_repro_mean, group = :food_level, lw = 2)
    display(plt_bestfit)
end

data.scalar["misc_traits"]
yhat_bestfit.scalar["misc_traits"]

[mode(p.untruncated) for p in priors.priors]
