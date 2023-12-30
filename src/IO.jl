function extract_colnames(c::R, k::Symbol) where R <: Real
    return k
end

function extract_colnames(c::AbstractVector, k::Symbol)
    return [Symbol("$(k)_$(i)") for i in 1:length(c)]
end

# TODO: this should not be needed any more - confirm that it can be removed
#function extract_colnames(c::AbstractMatrix, k::Symbol)
#    return vcat([[Symbol("$(k)_$(j)_$(i)") for i in 1:size(c)[1]] for j in 1:size(c)[2]]...)
#end

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

"""
Isolate the indicated PMoAs. 
I.e., turn off all PMoAs (including `h`!) except for those indicated in `pmoas`-Vector. <br>
This is done through the toxicokinetic rate constant.
$(TYPEDSIGNATURES)
"""
function isolate_pmoas!(deb::AbstractParams, pmoas::Vector{String})
    num_stressors = length(deb.k_D_G)
    deactivate = filter(x -> !(x in pmoas), ["G", "M", "A", "R"])
    for j in deactivate
        setfield!(deb, Symbol("k_D_$(j)"), zeros(num_stressors))
    end
    return nothing
end


function assert!(p::T) where T <: NamedTuple
    @assert length(p.deb.k_D_G) >= length(p.glb.C_W) "Length of k_D_G is not at least length of C_W"
    @assert length(p.deb.k_D_M) >= length(p.glb.C_W) "Length of k_D_M is not at least length of C_W"
    @assert length(p.deb.k_D_A) >= length(p.glb.C_W) "Length of k_D_A is not at least length of C_W"
    @assert length(p.deb.k_D_R) >= length(p.glb.C_W) "Length of k_D_R is not at least length of C_W"
    @assert length(p.deb.k_D_h) >= length(p.glb.C_W) "Length of k_D_h is not at least length of C_W"

    @assert length(p.deb.drc_functs_G) >= length(p.glb.C_W) "Length of drc_functs_G is not at least length of C_W"
    @assert length(p.deb.drc_functs_M) >= length(p.glb.C_W) "Length of drc_functs_G is not at least length of C_W"
    @assert length(p.deb.drc_functs_A) >= length(p.glb.C_W) "Length of drc_functs_G is not at least length of C_W"
    @assert length(p.deb.drc_functs_R) >= length(p.glb.C_W) "Length of drc_functs_G is not at least length of C_W"
    @assert length(p.deb.drc_functs_h) >= length(p.glb.C_W) "Length of drc_functs_G is not at least length of C_W"

    @assert length(p.deb.drc_params_G) >= length(p.glb.C_W) "Length of drc_params_G is not at least length of C_W"
    @assert length(p.deb.drc_params_M) >= length(p.glb.C_W) "Length of drc_params_M is not at least length of C_W"
    @assert length(p.deb.drc_params_A) >= length(p.glb.C_W) "Length of drc_params_A is not at least length of C_W"
    @assert length(p.deb.drc_params_R) >= length(p.glb.C_W) "Length of drc_params_R is not at least length of C_W"
    @assert length(p.deb.drc_params_h) >= length(p.glb.C_W) "Length of drc_params_h is not at least length of C_W"
end

