
#"""
#    function get_par_names(priors.HierchPriors)::Vector{Symbol}
#
#Get parameter names of a hierarhical prior structure in flattened form.
#"""
#function get_par_names(priors::HierchPriors)::Vector{Symbol}
#    
#    return vcat(
#        priors.hyperparam,
#        priors.groupparams,
#        priors.params
#    )
#end
#
#"""
#    linkfunction_zoomfactor(hypersample::Float64, numgroups::Int64)::Vector{Float64}
#
#Default link function for hierarchical model, assuming that the hyper parameter is a zoom factor.
#"""
#function linkfunction_zoomfactor(hypersample::Float64, numgroups::Int64)::Vector{Float64}
#    Zdist = Truncated(Normal(1, hypersample), 0, Inf)
#    return rand(Zdist, numgroups)
#end
#
#
#"""
#Structure for hierarchical priors.
#"""
#@with_kw mutable struct HierchPriors
#    hyperparam::Symbol
#    hyperprior::Distribution
#    groupparams::Vector{Symbol}
#    linkfunction::Function
#    params::Vector{Symbol}
#    priors::Vector{Distribution}
#
#    """
#        HierchPriors(
#            hyperprior::Pair{Symbol,D},
#            groupparams::Vector{Symbol}, 
#            priors::Vector{Pair{Symbol,D}} 
#            )
#
#    Iniitalize hierarhcical prior structure.
#
#    args:
#
#    - `hyperprior`: The name and distribution of the hyperprior as a Symbol/Distribution-pair
#    - `groupparams`: A list of names for the group parameters
#    - `priors`: The remaining parameters as a vector of Symbol/Distribution-pairs
#    """
#    function HierchPriors(
#        hyperprior::Pair{Symbol,D},
#        groupparams::Vector{Symbol}, 
#        priors::Vector{Pair{Symbol,D}};
#        linkfunction = linkfunction_zoomfactor
#        ) where D <: Distribution
#
#        hierchpriors = new()
#
#        hierchpriors.hyperparam = hyperprior.first
#        hierchpriors.hyperprior = hyperprior.second
#        hierchpriors.groupparams = groupparams
#        hierchpriors.params = [p.first for p in priors]
#        hierchpriors.priors = [p.second for p in priors]
#        hierchpriors.linkfunction = linkfunction
#        
#        return hierchpriors
#    end
#end

#"""
#    hierch_sample(priors::ABC.HierchPriors)::Vector{Float64}
#
#Take a sample from a hierarhcical prior structure, returning the sampled values 
#flattened as a `Vector{Float64}`.
#"""
#function hierch_sample(priors::ABC.HierchPriors)::Vector{Float64}
#    hypersample = rand(priors.hyperprior)
#    groupparam_samples = priors.linkfunction(hypersample, length(priors.groupparams))
#    param_samples = [rand(p) for p in priors.priors]
#
#    return vcat(hypersample, groupparam_samples, param_samples)
#end



"""
    initialize(
        n_pop::Int64, 
        defaultparams::AbstractParams, 
        priors::Tuple{Priors,DataFrame},
        simulator,
        distance,
        data::AbstractDataFrame
    )

Intialization of the SMC population for a hierarchical prior structure.
"""
#function initialize(
#    n_pop::Int64, 
#    defaultparams::Union{AbstractParams,AbstractParamCollection}, 
#    priors::HierchPriors,
#    simulator,
#    distance,
#    data::Any
#    )
#
#    particles = Vector{Vector{Float64}}(undef, n_pop) #[copy(defaultparams) for _ in 1:n_pop] # predefine vector of parameter samples
#    distances = Vector{Float64}(undef, n_pop)
#    weights = ones(n_pop) |> x -> x ./ sum(x)
#
#    @info("...Evaluating initial population...")
#
#    @threads for i in eachindex(particles)
#        particles[i] = hierch_sample(priors)
#        prediction = simulator(defaultparams, priors.params, particles[i])
#        distances[i] = distance(prediction, data)
#    end
#
#    return particles, distances, weights 
#
#end
