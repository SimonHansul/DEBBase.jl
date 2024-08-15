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
    function get_par_names(priors.HierchPriors)::Vector{Symbol}

Get parameter names of a hierarhical prior structure in flattened form.
"""
function get_par_names(priors::HierchPriors)::Vector{Symbol}
    
    return vcat(
        priors.hyperparam,
        priors.groupparams,
        priors.params
    )
end

"""
    linkfunction_zoomfactor(hypersample::Float64, numgroups::Int64)::Vector{Float64}

Default link function for hierarchical model, assuming that the hyper parameter is a zoom factor.
"""
function linkfunction_zoomfactor(hypersample::Float64, numgroups::Int64)::Vector{Float64}
    Zdist = Truncated(Normal(1, hypersample), 0, Inf)
    return rand(Zdist, numgroups)
end


"""
Structure for hierarchical priors.
"""
@with_kw mutable struct HierchPriors
    hyperparam::Symbol
    hyperprior::Distribution
    groupparams::Vector{Symbol}
    linkfunction::Function
    params::Vector{Symbol}
    priors::Vector{Distribution}

    """
        HierchPriors(
            hyperprior::Pair{Symbol,D},
            groupparams::Vector{Symbol}, 
            priors::Vector{Pair{Symbol,D}} 
            )

    Iniitalize hierarhcical prior structure.

    args:

    - `hyperprior`: The name and distribution of the hyperprior as a Symbol/Distribution-pair
    - `groupparams`: A list of names for the group parameters
    - `priors`: The remaining parameters as a vector of Symbol/Distribution-pairs
    """
    function HierchPriors(
        hyperprior::Pair{Symbol,D},
        groupparams::Vector{Symbol}, 
        priors::Vector{Pair{Symbol,D}};
        linkfunction = linkfunction_zoomfactor
        ) where D <: Distribution

        hierchpriors = new()

        hierchpriors.hyperparam = hyperprior.first
        hierchpriors.hyperprior = hyperprior.second
        hierchpriors.groupparams = groupparams
        hierchpriors.params = [p.first for p in priors]
        hierchpriors.priors = [p.second for p in priors]
        hierchpriors.linkfunction = linkfunction
        
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
