
# TODO: define a better SMC struct
#mutable struct SMC <: AbstractFitObject
#    data::Union{Dataset,} = "" # data can be reference to an object, function to generate quasi-observations or an object holding the data itsel
#    paramnames::Vector{Symbol} = [...] # parameter names
#    priors::Vector{Distribution} = [...] # parameter priors
#    hyperparam_names::Vector{Symbol} = [...] # hyperparameter names
#    hyperparam_priors::Vector{Distribution} = [...] # hyper priors
#    substructs::Union{Symbol,Vector{Symbol}} = [:pth, :pth, :spc, :glb,...] #
#    assignment::Vector{Function} = repeat([x -> x], length(paramnames)) # additional instructions for parameter assignments prior to simulation (e.g. back-transformation of log-transformed parameters)
#    paridx::OrderedDict{Symbol,Int64} = OrderedDict(zip(paramnames, eachindex(paramnames))) # parameter indexing to access priors etc. by name
#    epsilon_schedule::Vector{Float64} = Float64[] # threshold values
#    populations::Vector{DataFrame} = DataFrame[]
#    comptime::Union{Missing,Float64} = missing # computation time [s]
#    function SMC()
#    end
#end

@with_kw mutable struct SMCResult
    accepted::AbstractDataFrame = DataFrame()
    defaultparams::Any = []
    priors::Union{Priors,Tuple{Priors,DataFrame}} = Priors()
    simulator = nothing
    distance = nothing
    data::Any = nothing
    n_pop::Int64 = 0
    q_eps::Float64 = NaN
    k_max::Int64 = 0
    convergence_eps::Float64 = NaN
    time_of_execution::Union{String,DateTime} = ""
    comptime::Any = 0
    distance_schedule::Vector{Float64} = Float64[]
    converged::Bool = false # TODO: remove this in v1.0 (keeping it for now to maintain backwards compatability for v<1)
    errlog::DataFrame = DataFrame()
end


"""
initialize(
    n_pop::Int64, 
    defaultparams::Union{AbstractParams,AbstractParamCollection}, 
    priors::Priors,
    simulator,
    distance,
    data::AbstractDataFrame
    )

Initialization of a population if the priors are purely parameteric priors.
"""
function initialize_threaded(
    n_pop::Int64, 
    defaultparams::Union{AbstractParams,AbstractParamCollection}, 
    priors::Priors,
    simulator,
    distance,
    data::Any
    )
    particles = Vector{Vector{Float64}}(undef, n_pop) #[copy(defaultparams) for _ in 1:n_pop] # predefine vector of parameter samples
    distances = Vector{Float64}(undef, n_pop)
    weights = ones(n_pop) |> x -> x ./ sum(x)

    @info("...Evaluating initial population...")

    @threads for i in eachindex(particles)
        particles[i] = [rand(p) for p in priors.priors] 
        prediction = simulator(defaultparams, priors.params, particles[i])
        distances[i] = distance(prediction, data)
    end

    return particles, distances, weights 
end

"""
    initialize(
        n_pop::Int64, 
        defaultparams::AbstractParams, 
        priors::Tuple{Priors,DataFrame},
        simulator,
        distance,
        data::AbstractDataFrame
    )
    
Initialization of the SMC population if the given priors are a combination of parameteric priors and previously accepted particles. \\
Note that parameteric prior distributions should always be given for all parameters. 
For the parameters which also appear in the dataframe of accepted values
"""
function initialize_threaded(
    n_pop::Int64, 
    defaultparams::Union{AbstractParams,AbstractParamCollection}, 
    priors::Tuple{Priors,DataFrame},
    simulator,
    distance,
    data::Any
    )

    particles = Vector{Vector{Float64}}(undef, n_pop) # predefine vector of parameter samples
    distances = Vector{Float64}(undef, n_pop)
    weights = ones(n_pop) |> x -> x ./ sum(x)

    @info("...Evaluation of initial population...")

    @threads for i in 1:n_pop
        particles[i] = [rand(p) for p in priors[1].priors] # sample from parameteric priors
        
        acc_sample = posterior_sample(priors[2]) # get sample of accepted values and assign
        j = 0 # index counter for acc_sample
        for k in eachindex(particles[i]) # for all parameters
            if String(priors[1].params[k]) in DEBABC.get_par_names(priors[2]) # check if the value is supplied in the accepted values
                j += 1 # if so, increment counter
                particles[i][k] = acc_sample[j] # assign the value
            end
        end

        prediction = simulator(defaultparams, priors[1].params, particles[i])
        distances[i] = distance(prediction, data)
    end

    return particles, distances, weights 
end

function reject(
    particles::Vector{Vector{Float64}},
    distances::Vector{Float64}, 
    weights::Vector{Float64}, 
    q_eps::Float64
    ) 

    epsilon = quantile(distances[isfinite.(distances)], q_eps) # acceptance threshold
    mask = distances .< epsilon # indices of accepted samples
    accepted_particles = particles[mask] # accepted samples
    accepted_distances = distances[mask] # distances associated with accepted samples
    accepted_weights = weights[mask] # weights associated with accepted samples

    return accepted_particles, accepted_distances, accepted_weights, epsilon
end

"""
    gaussiankernel(old_value::Float64, scale::Float64, min::Float64, max::Float64)
Define a Gaussian perturbation kernel.
"""
function gaussiankernel(old_value::Float64, scale::Float64, min::Float64, max::Float64)
    return Truncated(Normal(old_value, scale), min, max)
end

"""
Calculate quantile-based parameter scale for scalar parameters.
"""
function calcscale(pvals::Vector{Float64}; qs::Tuple{Float64,Float64} = (.25, .75))
    return quantile(pvals, qs[2])  - quantile(pvals, qs[1])
end

"""
Calculate perturbation kernel scales from a Vector of accepted particles.
"""
function calculatescales(
    num_params::Int64, 
    parnames::Vector{Symbol}, 
    accepted_particles::Vector{Vector{Float64}}
    )
    scales = Vector(undef, num_params) # initialize scales
    for (k,parname) in enumerate(parnames)
        pvals = [a[k] for a in accepted_particles] #pvec(accepted_particles, parname) # extract parameter values
        scale = calcscale(pvals) # calculate kernel scale
        scales[k] = scale
    end
    return scales
end

#function getparam(particle::AbstractParams, paramname::Symbol)
#    if isvecparam(paramname, particle)
#        index = getparindex(paramname)
#        fieldname = getfieldname(paramname)
#        return getproperty(particle, fieldname)[index]
#    else
#        return getproperty(particle, paramname)
#    end
#end

"""
    resample(
        accepted_particles::Vector{Vector{Float64}}, 
        accepted_distances::Vector{Float64},
        accepted_weights::Vector{Float64},
        priors::Priors,
        n_pop::Int64,
        parnames::Vector{Symbol},
        num_params::Int64
        )

Resample from previously accepted particles, with account for SMC sampling weights.
"""
function resample(
    accepted_particles::Vector{Vector{Float64}}, 
    accepted_distances::Vector{Float64},
    accepted_weights::Vector{Float64},
    priors::Priors,
    n_pop::Int64,
    parnames::Vector{Symbol},
    num_params::Int64
    )
    
    scales = calculatescales(num_params, parnames, accepted_particles)
    idcs = sample(eachindex(accepted_particles), Weights(accepted_weights), n_pop) # indices of the re-sampled particles
    
    particles = [copy(accepted_particles[i]) for i in idcs] # resampled particles
    distances = [accepted_distances[i] for i in  idcs] # distances associated with resampled particles
    weights = [accepted_weights[i] for i in idcs] #accepted_weights[idcs] # weights associated with resampled particles

    for (i,particle) in enumerate(particles) # iterate over particles
        for (k,parname) in enumerate(parnames) # iterate over parameters
            old_value = particle[k] #getparam(particle, parname) # retrieve the old value
            kernel = gaussiankernel( # define the perturbation kernel
                    old_value, # make sure boundaries are respected
                    scales[k], 
                    minimum(priors.priors[k]), 
                    maximum(priors.priors[k])
                )
            new_value = rand(kernel) # sample the perturbed value
            particle[k] = new_value #assign!(particle, parname, new_value) # update value in the particle
            particles[i] = particle # update the particle
        end
    end

    return scales, idcs, particles, distances, weights
end

function resample(
    accepted_particles::Vector{Vector{Float64}}, 
    accepted_distances::Vector{Float64},
    accepted_weights::Vector{Float64},
    priors::Tuple{Priors,DataFrame},
    n_pop::Int64,
    parnames::Vector{Symbol},
    num_params::Int64
    )

    return resample(
        accepted_particles, 
        accepted_distances,
        accepted_weights,
        priors[1],
        n_pop,
        parnames,
        num_params
    )

end

"""
Calculate SMC sampling weights.
"""
function calculateweights_threaded(
    particles::Vector{Vector{Float64}}, 
    priors::Priors, 
    parnames::Vector{Symbol}, 
    accepted_particles::Vector{Vector{Float64}},
    accepted_idcs::Vector{Int64},
    accepted_weights::Vector{Float64},
    scales::AbstractVector
    )

    @info("...Calculating weights...")

    weights = Vector(undef, length(particles))
    
    @threads for i in eachindex(particles) # for every perturbed particle
        particle = particles[i] # get the particle
        priorprob = prod([pdf(prior, particle[k]) for (k,prior) in enumerate(priors.priors)]) # calculate prior probability
        denominator = 0. # start to calculate the denominator, which is a sum
        for (j,particle) in enumerate(particles) # for every particle
            kernelprob = 1 # start to calculate the probability of the particle according to the perturbation kernel
            old_particle = accepted_particles[accepted_idcs[j]] # get the unperturbed value
            old_weight = accepted_weights[accepted_idcs[j]] # get the old weight
            for (k,parname) in enumerate(parnames) # for every parameter in the particle
                old_value = old_particle[k] #getparam(old_particle, parname) # get the unperturbed value
                value = particle[k] #getparam(particle, parname) # get the perturbed value
                kernel = gaussiankernel( # define the perturbation kernel
                    old_value, # make sure boundaries are respected
                    scales[k], # use scales to define kernel
                    minimum(priors.priors[k]), # respect lower limit
                    maximum(priors.priors[k])  # respect upper limit
                )
                kernelprob *= pdf(kernel, value) # update the probability
            end
            denominator += old_weight * kernelprob # update the denominator
        end
        weights[i] = priorprob / denominator # update the weight
    end
    weights = weights ./ sum(weights) # normalize the weights

    return weights
end

function calculateweights_threaded(
    particles::Vector{Vector{Float64}}, 
    priors::Tuple{Priors,DataFrame}, 
    parnames::Vector{Symbol}, 
    accepted_particles::Vector{Vector{Float64}},
    accepted_idcs::Vector{Int64},
    accepted_weights::Vector{Float64},
    scales::AbstractVector
    )

    calculateweights_threaded(
        particles, priors[1], parnames, accepted_particles, accepted_idcs, accepted_weights, scales
    )

end

function evaluate_threaded(
    particles::Vector{Vector{Float64}},
    defaultparams::Union{AbstractParams,AbstractParamCollection},
    simulator,
    distance,
    priors::Priors,
    data::Any
    )

    @info("...Evaluating particles...")
    distances = Vector{Float64}(undef, length(particles))

    @threads for i in eachindex(particles)
        particle = particles[i]
        prediction = simulator(defaultparams, priors.params, particle)
        dist = distance(prediction, data)
        distances[i] = dist
    end

    return distances
end


function evaluate_threaded(
    particles::Vector{Vector{Float64}},
    defaultparams::Union{AbstractParams,AbstractParamCollection},
    simulator,
    distance,
    priors::Tuple{Priors,DataFrame},
    data::Any
    )

    return evaluate_threaded(
        particles,
        defaultparams,
        simulator,
        distance,
        priors[1],
        data
    )
end

#function check_convergence(k::Int64, distance_schedule::Vector{Float64}, convergence_eps::Float64)
#    distance_reldiff = (distance_schedule[k+1] - distance_schedule[k])/distance_schedule[k]
#    if distance_reldiff > 0
#        return false
#    elseif distance_reldiff >= convergence_eps
#        return true
#    else
#        return false
#    end
#end

get_par_names(priors::Priors) = priors.params # get parameter names from a parametric oject
get_par_names(priors::Tuple{Priors,DataFrame}) = unique(vcat([Symbol.(get_par_names(p)) for p in priors]...)) # get parameter names if a dataframe of accepted values is also supplied

@enum ParallelizationType sequential threaded distributed

"""
    SMC(
        priors::Union{Priors,Tuple{Priors,DataFrame}},
        defaultparams,
        simulator,
        distance,
        data::Any;
        n_pop::Int64 = 1000,
        q_eps::Float64 = 0.2,
        k_max::Int64 = 3,
        convergence_eps::Float64 = 0.1,
        paralellization_type = threaded,
        savetag = "smc",
        saveto::Union{String,False} = ""
        )

kwargs 

- `convergence_eps::Float64` : we consider SMC converged if the relative difference between two successive distance thresholds ϵ is less than `convergence_eps`
- `savedata::Bool` : whether to include the value of positional argument `data` in the output metadata. `false` by default
- `saveto::Bool` : where to store the `SMCResult` object on disc. `tempname()` by default.

"""
function SMC(
    priors::Union{Priors,Tuple{Priors,DataFrame}},
    defaultparams,
    simulator,
    distance,
    data::Any;
    n_pop::Int64 = 1000,
    q_eps::Float64 = 0.2,
    k_max::Int64 = 3,
    convergence_eps::Float64 = 0.1,
    paralellization_type = threaded,
    savetag = "smc",
    saveto::Union{String,Bool} = ""
    )
    

    @info("Executing SMC with $((k_max+1) * n_pop) samples on $(Threads.nthreads()) threads.")
    start = now() # record computation time

    let accepted_particles::Vector{Vector{Float64}},
        accepted_weights::Vector{Float64}, 
        k = 0, 
        converged = false,
        distance_schedule = Vector{Float64}(undef, k_max) # record the distance schedule
      
        param_names = get_par_names(priors)
        num_params = length(param_names)
        
        if paralellization_type == threaded
            particles, distances, weights = initialize_threaded(
                n_pop, defaultparams, priors, simulator, distance, data
                )
        else
            error("Paralellization types other than threaded yet to be implemented")
        end

        accepted_particles, accepted_distances, accepted_weights, epsilon = reject(
            particles, distances, weights, q_eps
            )

        push!(distance_schedule, epsilon)
        
        while (k < k_max) & (!converged)
            k += 1
            @info("#### Evaluating SMC step k = $(k) ####")

            scales, idcs, particles, distances, weights = resample(
                accepted_particles, accepted_distances, accepted_weights, priors, n_pop, param_names, num_params
                )

            if paralellization_type == threaded
                weights = calculateweights_threaded(
                    particles, priors, param_names, accepted_particles, idcs, accepted_weights, scales
                    )
            else
                error("Paralellization types other than threaded yet to be implemented")
            end
    

            if paralellization_type == threaded
                distances = evaluate_threaded(
                    particles, defaultparams, simulator, distance, priors, data
                    )
            else
                error("Paralellization types other than threaded yet to be implemented")
            end
        
            accepted_particles, accepted_distances, accepted_weights, epsilon = reject(
                particles, distances, weights, q_eps
                )
            
            distance_schedule[k] = epsilon
            #converged = check_convergence(k, distance_schedule, convergence_eps)
        end

        @info("SMC reached k_max")
       
        # assemble dataframe of accepted particles

        accepted = [[particle[k] for k in eachindex(get_par_names(priors))] for particle in accepted_particles] |> # convert Vector of particles to nested Vector of parameter values
        x -> hcat(x...)' |> # convert nested Vector to Matrix
        x -> DataFrame(x, param_names) # convert Matrix to DataFrame
        accepted[!,:distance] .= accepted_distances # record associated distances
        accepted[!,:weight] .= accepted_weights # record SMC weights

        # assemble metadata

        fit = SMCResult(
            accepted = accepted,
            defaultparams = defaultparams,
            priors = priors,
            simulator = simulator,
            distance = distance,
            data = data,
            n_pop = n_pop,
            q_eps = q_eps,
            k_max = k_max,
            convergence_eps = convergence_eps,
            time_of_execution = start,
            comptime = now() - start, 
            distance_schedule = distance_schedule,
            converged = converged
            )

        if saveto != false
            CSV.write(joinpath("$saveto", "$(savetag)_accepted.csv"), accepted)
            CSV.write(joinpath("$saveto", "$(savetag)_info.csv"), DataFrame(
                defaultparams = defaltaparams,
                priors = priors,
                n_pop = n_pop,
                q_eps = q_eps,
                k_max = k_max,
                time_of_execution = start,
                comptime = fit.comptime

            ))
        end

        return fit
    end
end
