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


