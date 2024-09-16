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
