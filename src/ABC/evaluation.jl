"""
Weighted Kernel Density Estimate. 
"""
function wkde(X::Vector{R}, ω = Vector{R}) where R <: Real
    U = kde(X, weights = Weights(ω))
    return U
end

function wkde(accepted::DataFrame, param::Symbol)
    return wkde(accepted[:,param], accepted[:,:weight])
end

function wkde(X::Matrix{R}, ω = Vector{R}) where R <: Real
    U = kde(X, weights = Weights(ω))
    return U 
end

function wkde(accepted, param1::Symbol, param2::Symbol)
    X = Matrix(hcat([accepted[:,param1], accepted[:,param2]]...))
    return wkde(X, accepted[:,:weight])
end


function get_par_names(accepted)
    return names(accepted)[1:length(names(accepted))-2]
end

"""
    summarize_accepted(
        accepted::D; 
        ffun = fround
        ) where D <: AbstractDataFrame  
Compute summary statistics from a data frame of accepted particles. 
Values are formatted using ``ffun``. \n
Use ``ffun = x -> x`` to suppress any kind of formatting.
"""
function summarize_accepted(
    accepted::D; 
    ffun = fround
    ) where D <: AbstractDataFrame  

    ptest = accepted[accepted.distance.==minimum(accepted.distance),:]
    parnames = get_par_names(accepted) |> x -> Symbol.(x)

    best_fits = []
    means = []
    medians = []
    q05 = []
    q25 = []
    q75 = []
    q95 = []

    for par in parnames
        let x = accepted[:,par]
            push!(best_fits, ptest[1,par])
            push!(means, mean(x))
            push!(medians, median(x))
            push!(q05, quantile(x, 0.05))
            push!(q25, quantile(x, 0.25))
            push!(q75, quantile(x, 0.75))
            push!(q95, quantile(x, 0.95))
        end

    end

    summary = DataFrame(
        param = parnames, 
        bestfit = ffun.(best_fits), 
        mean = ffun.(means),
        median = ffun.(medians),
        q05 = ffun.(q05), 
        q25 = ffun.(q25), 
        q75 = ffun.(q75),
        q95 = ffun.(q95)
    )

    return summary
end

"""
    ppc(defaultparams::AbstractParams, simulator, accepted::AbstractDataFrame, priors::Priors; n_samples = 1000) 

Run a posterior predictive check, given 

- `defaultparam::AbstractParams`: default parameters used in Simulator
- `simulator`: simulator function
- `accepted::AbstractDataFrame`> accepted particles returned by SMC

This function assumes that `simulator` returns a `DataFrame`.

kwargs

- `n_samples`: number of samples from `accepted_particles` to evaluate
"""
function ppc(defaultparams::Union{AbstractParams,AbstractParamCollection}, simulator, accepted::AbstractDataFrame, priors::Priors; n_samples = 1000) 
    @info("Running posterior predictive check with $n_samples samples on $(Threads.nthreads()) threads")
    yhat = Vector{DataFrame}(undef, n_samples) # predictions

    @threads for i in 1:n_samples # for the given number of samples
        sample = posterior_sample(accepted)
        yhat_i = simulator(defaultparams, priors.params, sample) # evaluate the sample 
        yhat_i[!,:n_sample] .= i # add ID column
        yhat[i] = yhat_i # collect results
    end

    yhat = vcat(yhat...) # convert Vector of DataFrames to DataFrame

    return yhat
end

"""
    ppc(res::SMCResult; n_samples = 1000, kwargs...)
    
Run a posterior predictive check, given SMC results `res`.
"""
function ppc(res::SMCResult; n_samples = 1000, kwargs...)
    return ppc(res.defaultparams, res.simulator, res.accepted, res.priors; n_samples = n_samples, kwargs...)
end



