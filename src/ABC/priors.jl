mutable struct Priors
    params::Vector{Symbol}
    priors::Vector{Distribution}

    """
        Priors(args::Pair...)
        
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


    """
        Priors(params, priors)

    Initialize priors from a vector of parameter names and prior distributions, respectively.
    """
    function Priors(params, priors)
        return new(params, priors)
    end
end


