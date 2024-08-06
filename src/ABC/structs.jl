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

function get_par_names(numgroups::Int64, groupnames::Nothing, hyperpriors::Vector{Pair{Symbol,D}}, priors::Vector{Pair{Symbol,D}}) where D <: Distribution

end

@with_kw mutable struct HierchPriors
    hyperparams::Vector{Symbol}
    hyperpriors::Vector{Distribution}
    numgroups::Int64
    groupnames::Union{Nothing,Vector{Symbol}} = nothing
    params::Vector{Symbol}
    priors::Vector{Distribution}

    function HierchPriors(
        numgroups, groupnames; 
        hyperpriors::Vector{Pair{Symbol,D}}, 
        priors::Vector{Pair{Symbol,D}} 
        )

        parnames = get_par_names(numgroups, groupnames, hyperpriors, priors)


        
    end
end



@with_kw mutable struct SMCResult
    accepted::AbstractDataFrame
    defaultparams::AbstractParams
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
