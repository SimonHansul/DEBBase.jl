"""
    wkde(X::Vector{R}, ω = Vector{R}) where R <: Real

Compute a weighted Kernel Density Estimate of values in `X` with weights `ω`. 
"""
function wkde(X::Vector{R}, ω = Vector{R}) where R <: Real
    U = kde(X, weights = Weights(ω))
    return U
end

"""
    function wkde(accepted::DataFrame, param::Symbol)

Compute a weighted kernel density estimate for parameter `param` in accepted particles `accepted`.
"""
function wkde(accepted::DataFrame, param::Symbol)
    return wkde(accepted[:,param], accepted[:,:weight])
end

"""
    function wkde(X::Matrix{R}, ω = Vector{R}) where R <: Real

Compute a weighted kernel density estimate on Matrix `X`.
"""
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
Values are formatted using function ``ffun``. \n
Use ``ffun = x -> x`` to suppress any kind of formatting. 
By default, values are rounded to 2 significant digits and formatted as strings using `fround`.
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
    ppc(defaultparams::Any, simulator, accepted::AbstractDataFrame, priors::Priors; n_samples = 1000) 

Compute posterior predictions for posterior predictive check.

- `defaultparam::Any`: default parameters used in simulator
- `simulator`: simulator function with signature `simulator(defaultparams, priors.params, sample)`, where 
- `accepted::AbstractDataFrame`: accepted particles as returned by SMC

This function assumes that `simulator` returns a `DataFrame`.

kwargs

- `n_samples = 1000`: number of samples from `accepted_particles` to evaluate
"""
function ppc(defaultparams::Any, simulator, accepted::AbstractDataFrame, priors::Priors; n_samples = 1000) 
    @info("Running posterior predictive check with $n_samples samples on $(Threads.nthreads()) threads")
    sim = Vector{DataFrame}(undef, n_samples) # predictions

    @threads for i in 1:n_samples # for the given number of samples
        sample = posterior_sample(accepted)
        sim_i = simulator(defaultparams, priors.params, sample) # evaluate the sample 
        sim_i[!,:n_sample] .= i # add ID column
        sim[i] = sim_i # collect results
    end

    sim = vcat(sim...) # convert Vector of DataFrames to DataFrame

    return sim
end

"""
    ppc(fit::SMCResult; n_samples = 1000, kwargs...)
    
Run a posterior predictive check, given SMCResult object `fit`.
"""
function ppc(fit::SMCResult; n_samples = 1000, kwargs...)
    return ppc(fit.defaultparams, fit.simulator, fit.accepted, fit.priors; n_samples = n_samples, kwargs...)
end



