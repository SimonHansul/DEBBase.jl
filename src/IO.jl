function extract_colnames(c::R, k::Symbol) where R <: Real
    return k
end

function extract_colnames(c::AbstractVector, k::Symbol)
    return [Symbol("$(k)_$(i)") for i in 1:length(c)]
end

function extract_colnames(c::AbstractMatrix, k::Symbol)
    return vcat([[Symbol("$(k)_$(j)_$(i)") for i in 1:size(c)[1]] for j in 1:size(c)[2]]...)
end

"""
Extract the column names of output data frame from a ODE solution object. 
This function assumes that ComponentArrays is used.
"""
function extract_colnames(sol::O)::Vector{Symbol} where {O <:ODESolution}
    let u = sol.u[1]
        colnames = []

        for k in keys(u)
            push!(colnames, extract_colnames(u[k], k))
        end

        colnames = vcat([:t], colnames...)
        return colnames
    end
end

"""
Convert ODE solution object to output data frame.
$(TYPEDSIGNATURES)
"""
function sol_to_df(sol::O)::DataFrame where O <: ODESolution
    simout = DataFrame(
        hcat(sol.t, hcat(sol.u...)'), # ODE output converted to wide matrix
        extract_colnames(sol)
    )
    return simout
end