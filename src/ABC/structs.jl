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
