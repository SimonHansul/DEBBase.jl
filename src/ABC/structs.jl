mutable struct Priors
    params::Vector{Symbol}
    priors::Vector{Distribution}

    """
    Initialize priors with a sequence of Symobl/Distribution pairs.
    """
    function Priors(args::Pair...)
        params = Symbol[]
        priors = Distribution[]

        for pair in args
            push!(params, pair.first)
            push!(priors, pair.second)
        end

        return new(params, priors)
    end

    function Priors(arg::Vector{Pair{Symbol,D}}) where D <: Distribution
        params = Symbol[]
        priors = Distribution[]

        for pair in arg
            push!(params, pair.first)
            push!(priors, pair.second)
        end

        return new(params, priors)
    end
end

"""
get_par_names(
    numgroups::Int64, 
    groupnames::Nothing, 
    hyperpriors::Vector{Pair{Symbol,D}}, 
    priors::Vector{Pair{Symbol,D}}
    ) where D <: Distribution

Infer parameter names of a hierarchical model, where no group names were given. 
Groups will be assigned numbers.
"""
function get_par_names(
    numgroups::Int64, 
    groupnames::Nothing, 
    hyperpriors::Vector{Pair{Symbol,D}}, 
    priors::Vector{Pair{Symbol,D}}
    ) where D <: Distribution

    # "normal" parameter names
    parnames = [p.first for p in priors]
    groupnames = string.(1:numgroups)

    # parameter for each group
    for (suffix,hyperprior) in zip(groupnames,hyperpriors)
        parname = String(hyperprior.first)*"_"*suffix
        pushfirst!(parnames, parname)
    end

    return parnames
end

"""
    linkfunction_zoomfactor(hypersamples::Vector{Float64}, numgrups::Int64)

Default link function for hierarchical model, assuming that the hyper parameter is a zoom factor.
"""
function linkfunction_zoomfactor(hypersamples::Vector{Float64}, numgroups::Int64)
    Zdist = rand(Truncated(Normal(1, hypersamples[1]), 0, Inf))
    return rand(Zdist, numgroups)
end

"""
Structure for hierarchical priors.
"""
@with_kw mutable struct HierchPriors
    hyperparams::Vector{Symbol}
    hyperpriors::Vector{Distribution}
    groupparams::Vector{Symbol}
    linkfunction::Function
    params::Vector{Symbol}
    priors::Vector{Distribution}

    function HierchPriors(
        hyperpriors::Vector{Pair{Symbol,D}},
        groupparams::Vector{Symbol}, 
        priors::Vector{Pair{Symbol,D}} 
        )

        hierchpriors = new()

        hierchpriors.hyperparams = [p.first for p in hyperpriors]
        hierchpriors.hyperpriors = [p.second for p in hyperpriors]
        hierchpriors.groupparams = groupparams
        hierchpriors.params = [p.first for p in priors]
        hierchpriors.priors = [p.second for p in priors]
        
        return hierchpriors
    end
end



@with_kw mutable struct SMCResult
    accepted::AbstractDataFrame
    defaultparams::Any
    priors::Union{Priors,Tuple{Priors,DataFrame}}
    simulator
    distance
    data::Any
    n_pop::Int64
    q_eps::Float64
    k_max::Int64
    convergence_eps::Float64
    time_of_execution::DateTime
    comptime
    distance_schedule::Vector{Float64}
    converged::Bool
end
