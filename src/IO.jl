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
Trim TKTD parameters to the smallest number of stressors indicated for any PMoA. <br>
For example `k_D_G = [1., 0.], k_D_M = [1.]` will be trimmed to `k_D_G = [1.], k_D_M = [1.]`.
$(TYPEDSIGNATURES)
"""
function trim!(deb::AbstractParams)
    num_stressors = min(
        length(deb.k_D_G),
        length(deb.k_D_M),
        length(deb.k_D_A),
        length(deb.k_D_R),
        length(deb.k_D_h)
    )

    deb.k_D_G = deb.k_D_G[1:num_stressors]
    deb.k_D_M = deb.k_D_M[1:num_stressors]
    deb.k_D_A = deb.k_D_A[1:num_stressors]
    deb.k_D_R = deb.k_D_R[1:num_stressors]
    deb.k_D_h = deb.k_D_h[1:num_stressors]

    deb.drc_functs_G = deb.drc_functs_G[1:num_stressors]
    deb.drc_functs_M = deb.drc_functs_M[1:num_stressors]
    deb.drc_functs_A = deb.drc_functs_A[1:num_stressors]
    deb.drc_functs_R = deb.drc_functs_R[1:num_stressors]
    deb.drc_functs_h = deb.drc_functs_h[1:num_stressors]

    deb.drc_params_G = deb.drc_params_G[1:num_stressors]
    deb.drc_params_M = deb.drc_params_M[1:num_stressors]
    deb.drc_params_A = deb.drc_params_A[1:num_stressors]
    deb.drc_params_R = deb.drc_params_R[1:num_stressors]
    deb.drc_params_h = deb.drc_params_h[1:num_stressors]
end


"""
Isolate the indicated PMoAs for chemical stressor `z`. 
That means, turn off all PMoAs (including lethal effects `h`) except for those indicated in `pmoas`-Vector. 
This is done through the toxicokinetic rate constant.
$(TYPEDSIGNATURES)
"""
function isolate_pmoas(deb::AbstractParams, pmoas::Vector{String}, z::Int64)::AbstractParams
    deactivate = filter(x -> !(x in pmoas), ["G", "M", "A", "R", "h"])
    for j in deactivate
        let fieldname = Symbol("k_D_$(j)")
            k_D = getfield(deb, fieldname)
            k_D[z] = 0.
            setfield!(deb, fieldname, k_D)
        end
    end
    return deb
end

function isolate_pmoas(deb::AbstractParams, pmoas::Vector{String})::AbstractParams
    deactivate = filter(x -> !(x in pmoas), ["G", "M", "A", "R", "h"])
    for j in deactivate
        let fieldname = Symbol("k_D_$(j)")
            k_D = getfield(deb, fieldname)
            k_D .= 0.
            setfield!(deb, fieldname, k_D)
        end
    end
    return deb
end

function isolate_pmoas!(deb::AbstractParams, pmoas::Vector{String}, z::Int64)::Nothing
    deb = isolate_pmoas(deb, pmoasm, z)
    return nothing
end

function isolate_pmoas!(deb::AbstractParams, pmoas::Vector{String})::Nothing
    deb = isolate_pmoas(deb, pmoas)
    return nothing
end

"""
Raise assertion errors
$(TYPEDSIGNATURES)
"""
function assert!(p::AbstractParamCollection)
    @assert length(p.deb.k_D_G) >= length(p.glb.C_W) "Length of k_D_G is not at least length of C_W"
    @assert length(p.deb.k_D_M) >= length(p.glb.C_W) "Length of k_D_M is not at least length of C_W"
    @assert length(p.deb.k_D_A) >= length(p.glb.C_W) "Length of k_D_A is not at least length of C_W"
    @assert length(p.deb.k_D_R) >= length(p.glb.C_W) "Length of k_D_R is not at least length of C_W"
    @assert length(p.deb.k_D_h) >= length(p.glb.C_W) "Length of k_D_h is not at least length of C_W"

    @assert length(p.deb.drc_functs_G) >= length(p.glb.C_W) "Length of drc_functs_G is not at least length of C_W"
    @assert length(p.deb.drc_functs_M) >= length(p.glb.C_W) "Length of drc_functs_M is not at least length of C_W"
    @assert length(p.deb.drc_functs_A) >= length(p.glb.C_W) "Length of drc_functs_A is not at least length of C_W"
    @assert length(p.deb.drc_functs_R) >= length(p.glb.C_W) "Length of drc_functs_R is not at least length of C_W"
    @assert length(p.deb.drc_functs_h) >= length(p.glb.C_W) "Length of drc_functs_h is not at least length of C_W"

    @assert length(p.deb.drc_params_G) >= length(p.glb.C_W) "Length of drc_params_G is not at least length of C_W"
    @assert length(p.deb.drc_params_M) >= length(p.glb.C_W) "Length of drc_params_M is not at least length of C_W"
    @assert length(p.deb.drc_params_A) >= length(p.glb.C_W) "Length of drc_params_A is not at least length of C_W"
    @assert length(p.deb.drc_params_R) >= length(p.glb.C_W) "Length of drc_params_R is not at least length of C_W"
    @assert length(p.deb.drc_params_h) >= length(p.glb.C_W) "Length of drc_params_h is not at least length of C_W"
end

"""
Calculate relative response, given scalar values.
"""
function rr(x::R, x_ref::R)::Float64 where R <: Real
    if x_ref == 0.
        return 1.
    else
        return x / x_ref
    end
end

function rr(x::R, x_ref::Missing)::Missing where R <: Real
    return missing
end


"""
Calculate the relative responses. \n
Positional arguments: 
- `res::AbstractDataFrame`: results
- `response_vars::Vector{Symbol}`: response variables for which to calculate the relative responses
- `treatment_var::Symbol`: Column indicating the treatment. Column values can be numerical or categorical, but `identify_control` kwarg has to be specified in the latter case

Keyword arguments: 
- `groupby_vars::Vector{Symbol}`: relative response will be conditioned on these variables (e.g. time, separate experiments...). Empty by default.
- `identify_control`: function used to identify reference values from `treatment_var`. By default, this is `minimum()` (assuming that column values in `treatment_var` are numerical).
"""
function relative_response(
    res::D, 
    response_vars::Vector{Symbol},
    treatment_var::Symbol; 
    groupby_vars::Vector{Symbol} = Symbol[],
    identify_control = minimum
    ) where D <: AbstractDataFrame

    #=
    Calculation of the conditional control mean
    =#
    reference = res[res[:,treatment_var] .== identify_control(res[:,treatment_var]),:] |> # extract the control
    x -> groupby(x, groupby_vars) |> # group by conditioning variables
    x -> combine(x) do df 
        refvals = DataFrame() # (sub-)dataframe of reference values
        for var in response_vars # iterate over response values
            var_ref = Symbol(String(var) * "_ref") # get the reference column name
            refvals[!,var_ref] = [mean(df[:,var])] # calculate the conditional control mean
        end
        return refvals
    end

    #=
    Calculation of the relative response
    =#
    res = leftjoin(res, reference, on = groupby_vars) # add references as new columns
    for var in response_vars # for every response variable
        y_var = Symbol("y_" * String(var)) # get the relative response column name
        var_ref = Symbol(String(var) * "_ref") # get the reference column name
        res[!,y_var] = [rr(row[var], row[var_ref]) for row in eachrow(res)] # calculate the relative response
        select!(res, Not(var_ref)) # drop the reference column
    end
    return res
end



"""
Set parameters with common prefix to the value of a reference parameter. 
E.g. 

```
set_equal!(deb, :Idot_max_rel, :lrv)
```

sets all parameters whose names start with `Idot_max_rel` equal to `Idot_max_rel_lrv`.
"""
function set_equal!(deb::AbstractParams, prefix::SS, ref_suffix::SS) where SS <: Union{Symbol,String}
    prefix = String(prefix)
    ref_suffix = String(ref_suffix)

    ref_param = prefix*"_"*ref_suffix
    paramnames = String[String.(fieldnames(typeof(deb)))...]
    filter!(x -> occursin(String(prefix), x), paramnames)
    filter!(x -> x != String(ref_param), paramnames)

    for param in paramnames
        setproperty!(deb, Symbol(param), Symbol(ref_param))
    end
end
