


"""
Sample from priors and assign to existing instance of parameter structure.
"""
function rand!(p::AbstractParams, priors::Priors)
    for (param,prior) in zip(priors)
        value = rand(prior)
        assign!(p, param, value)
    end
end

"""
    rand(defaultparams::DataType, priors::Priors)

Sample from priors and assign to new instance of default parameter structure.
"""
function rand(defaultparams::DataType, priors::Priors)
    newparams = defaultparams()
    rand!(newparams, priors)
    return newparams
end

function rand(defaultparams::A, priors::Priors) where A <: AbstractParams
    newparams = copy(defaultparams)
    rand!(newparams, priors)
    return newparams
end

function rand(priors::Priors)
    return [rand(p) for p in priors.priors]
end


"""
    posterior_sample(accepted::DataFrame; reserved_colnames::Vector{String} = ["distance", "weight", "model", "chain"])

Take posterior sample from a data frame of accepted values.
"""
function posterior_sample(accepted::DataFrame; reserved_colnames::Vector{String} = ["distance", "weight", "model", "chain"])::Vector{Float64}
    ω =  accepted.weight
    selectcols = filter(x -> !(x in reserved_colnames), names(accepted)) 
    sampled_values = accepted[sample(axes(accepted, 1), Weights(ω)),selectcols]
    return Vector{Float64}(sampled_values)
end

posterior_sample(res::SMCResult; kwargs...) = posterior_sample(res.accepted; kwargs...)


"""
Take posterior sample from a data frame of accepted values and assign to parameter structure. 
ID cols are at the end of the data frame and will be omitted from the sample, according to `num_idcols`.
"""
function posterior_sample!(p::AbstractParams, accepted::DataFrame; reserved_colnames::Vector{String} = ["distance", "weight", "model", "chain"])
    ω =  accepted.weight
    selectcols = filter(x -> !(x in reserved_colnames), names(accepted)) #1:length(names(accepted)) - num_idcols
    sampled_values = accepted[sample(axes(accepted, 1), Weights(ω)),selectcols]
    param_names = names(sampled_values)
    assign!(p, param_names, sampled_values)
end


#"""
#    posterior_sample(defparams::AbstractParams, accepted::AbstractDataFrame)
#Take posterior sample from `accepted` and assign to a copy of `defparams`.
#"""
#function posterior_sample(defparams::AbstractParams, accepted::AbstractDataFrame)
#    theta = copy(defparams)
#    DEBABC.posterior_sample!(theta, accepted)
#    return theta
#end


"""
    bestfit(defparams::AbstractParams, accepted::AbstractDataFrame)
Get the best fit from `accepted` (particle with minimum distance) and assign to a copy of `defparams`.
"""
function bestfit(defparams::AbstractParams, accepted::AbstractDataFrame)
    return posterior_sample(accepted[accepted.distance.==minimum(accepted.distance),:])
end
