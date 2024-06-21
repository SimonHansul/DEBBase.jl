function extract_colnames(c::R, k::Symbol) where R <: Real
    return k
end

function extract_colnames(c::AbstractVector, k::Symbol)
    return [Symbol("$(k)_$(i)") for i in 1:length(c)]
end

"""
    extract_colnames(u::ComponentVector)
Extract the column names of an output dataframe from an ODE solution object.
"""
function extract_colnames(u::ComponentVector)
    colnames = []
    for k in keys(u)
        push!(colnames, extract_colnames(u[k], k))
    end
    return vcat(colnames...)
end

"""
    extract_colnames(sol::O)::Vector{Symbol} where {O <:ODESolution}

Extract the column names of an output dataframe from an ODE solution object.
"""
extract_colnames(sol::ODESolution)::Vector{Symbol} = vcat([:t], extract_colnames(sol.u[1]))


