import Base:zip
"""
    zip(priors::Priors)
Simultaneously iterate over parameter names and prior distributions.
"""
function zip(priors::Priors)
    Base.zip(priors.params, priors.priors)
end

"""
Given parameter `param`, convert vector of parameter structures containg `param` into vector of values of `param`.
"""
function pvec(particles::Vector{P}, param::Symbol) where P <: AbstractParams

    paramvec = Vector{Float64}(undef, length(particles)) # get a flattened vector of parameter values)

    for (i,val) in enumerate(particles)
        if isvecparam(param, particles[1]) # account for parameters which are vector elements
            index = getparindex(param) # get index of the element
            fieldname = getfieldname(param) # get the valid param struct fieldname
            paramvec[i] = getproperty(val, fieldname)[index] # update vector index
        else # for parameters which are scalars, 
            paramvec[i] = getproperty(val, param) # simply update the field
        end
    end
    return paramvec
end

import Base:get
"""
    get(priors::Priors, param::Symbol)

Get prior distribution from priors object. 

## Examples:
    prior_kM = get(priors, :k_M)

"""
function get(priors::Priors, paramname::Symbol)
    return priors.priors[priors.params .== paramname][1]
end

"""
    isvecparam(paramname::Symbol, particle::AbstractParams)::Bool

Check whether parameter represents a Vector element. 
We assume that the indices of such parameters are indicated by a suffix following the form 
`param_i`. E.g. `k_D_G_1` for the first value in the vector of `k_D_G`-values. 
`isvecparam` checks whether such a field occurs in `AbstractParams` and whether it is a subtype of `AbstractVector`.
"""
function isvecparam(paramname::Symbol, particle::AbstractParams)::Bool 
    paramstem = join(split(String(paramname), "_")[1:end-1], "_") 
    
    # first check whether the stem of the parameter occurs in the fieldnames of params, 
    # then check whether it is a Vector
    return (Symbol(paramstem) in fieldnames(typeof(particle))) && (typeof(getproperty(particle, Symbol(paramstem))) <: AbstractVector)
end

"""
    getparindex(paramname::Symbol)::Int64

Infer index from a parameter name.

## Examples 

    getparindex(:k_D_G_1) # returns 1, because parameter name ends on _1
    getparindex(:k_D_G) # will fail, because no index indicated in suffix
"""
getparindex(paramname::Symbol)::Int64 = parse(Int64, split(String(paramname), "_")[end])

"""
    getfieldname(paramname::Symbol)::Symbol
Extract a valid param struct fieldname from a parameter name by dropping the suffix. 

## Examples
    getfieldname(:k_D_A_1) # returns :k_D_A
"""
getfieldname(paramname::Symbol)::Symbol = join(split(String(paramname), "_")[1:end-1], "_") |> Symbol

"""
    assign!(particle::AbstractParams, paramname::Symbol, value::Float64; assignment_instructions::Nothing)

Assign a sample to a param struct. 
This accounts for the possibility that some parameters might be stored in vectors (e.g. TKTD parameters where each vector element corresponds to a chemical). 
If kwarg `assignment_instructions` is `Nothing`, we can ignore it.
"""
function assign!(particle::AbstractParams, paramname::Symbol, value::Float64; assignment_instructions::Nothing)
    if isvecparam(paramname, particle) # does the parameter name indicate an index?
        index = getparindex(paramname) # extract index as integer
        fieldname = getfieldname(paramname) # extract the param struct fieldname
        vectorparam = copy(getproperty(particle, fieldname)) # get the vector parameter
        setindex!(vectorparam, value, index)
        setproperty!(particle, fieldname, vectorparam) # update the vector at the given index
    else # otherwise, 
        setproperty!(particle, paramname, value) # simply update the property
    end
end

"""
assign!(particle::AbstractParams, paramname::Symbol, value::Float64; assignment_instructions::Dict)

Assign a sample to a param struct. 
This accounts for the possibility that some parameters might be stored in vectors (e.g. TKTD parameters where each vector element corresponds to a chemical). 
If kwarg `assignment_instructions` is `Nothing`, we can ignore it.
"""
function assign!(particle::AbstractParams, paramname::Symbol, value::Float64; assignment_instructions::Dict)
    if isvecparam(paramname, particle) # does the parameter name indicate an index?
        index = getparindex(paramname) # extract index as integer
        fieldname = getfieldname(paramname) # extract the param struct fieldname
        vectorparam = copy(getproperty(particle, fieldname)) # get the vector parameter
        setindex!(vectorparam, value, index) # update the vector at the given index 
        setproperty!(particle, fieldname, vectorparam) # the property of the param struct, i.e. the whole vector
    else # otherwise, 
        if paramname in keys(assignment_instructions) # if an assignment instruction is given for this parameter, 
            assignment_instructions[paramname](particle, value) # apply it
        else # otherwise
            setproperty!(particle, paramname, value) # simply update the property
        end
    end
end

"""
    assign!(particle::AbstractParams, paramnames::Union{Vector{Symbol},Vector{String}}, pvals::Union{AbstractVector,DataFrameRow})

Assign new values to a parameter structure, given a list of names and values.
"""
function assign!(particle::AbstractParams, paramnames::Union{Vector{Symbol},Vector{String}}, pvals::Union{AbstractVector,DataFrameRow})
    for (param,value) in zip(paramnames,pvals)
        param = Symbol(param)
        assign!(particle, param, value)
    end
end