function extract_colnames(c::R, k::Symbol) where R <: Real
    return k
end

function extract_colnames(c::AbstractVector, k::Symbol)
    return [Symbol("$(k)_$(i)") for i in 1:length(c)]
end

# TODO: this should not be needed any more - confirm that it can be removed
#function extract_colnames(c::AbstractMatrix, k::Symbol)
#    return vcat([[Symbol("$(k)_$(j)_$(i)") for i in 1:size(c)[1]] for j in 1:size(c)[2]]...)
#endc:\Users\simon\Downloads\CV_2024_03.pd

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
function trim!(spc::AbstractParams)
    num_stressors = min(
        length(spc.k_D_G),
        length(spc.k_D_M),
        length(spc.k_D_A),
        length(spc.k_D_R),
        length(spc.k_D_h)
    )

    spc.k_D_G = spc.k_D_G[1:num_stressors]
    spc.k_D_M = spc.k_D_M[1:num_stressors]
    spc.k_D_A = spc.k_D_A[1:num_stressors]
    spc.k_D_R = spc.k_D_R[1:num_stressors]
    spc.k_D_h = spc.k_D_h[1:num_stressors]

    spc.drc_functs_G = spc.drc_functs_G[1:num_stressors]
    spc.drc_functs_M = spc.drc_functs_M[1:num_stressors]
    spc.drc_functs_A = spc.drc_functs_A[1:num_stressors]
    spc.drc_functs_R = spc.drc_functs_R[1:num_stressors]
    spc.drc_functs_h = spc.drc_functs_h[1:num_stressors]

    spc.drc_params_G = spc.drc_params_G[1:num_stressors]
    spc.drc_params_M = spc.drc_params_M[1:num_stressors]
    spc.drc_params_A = spc.drc_params_A[1:num_stressors]
    spc.drc_params_R = spc.drc_params_R[1:num_stressors]
    spc.drc_params_h = spc.drc_params_h[1:num_stressors]
end


"""
Isolate the indicated PMoAs for chemical stressor `z`. 
That means, turn off all PMoAs (including lethal effects `h`) except for those indicated in `pmoas`-Vector. 
This is done through the toxicokinetic rate constant.
$(TYPEDSIGNATURES)
"""
function isolate_pmoas(spc::AbstractParams, pmoas::Vector{String}, z::Int64)::AbstractParams
    deactivate = filter(x -> !(x in pmoas), ["G", "M", "A", "R", "h"])
    for j in deactivate
        let fieldname = Symbol("k_D_$(j)")
            k_D = getfield(spc, fieldname)
            k_D[z] = 0.
            setfield!(spc, fieldname, k_D)
        end
    end
    return spc
end

function isolate_pmoas(spc::AbstractParams, pmoas::Vector{String})::AbstractParams
    deactivate = filter(x -> !(x in pmoas), ["G", "M", "A", "R", "h"])
    for j in deactivate
        let fieldname = Symbol("k_D_$(j)")
            k_D = getfield(spc, fieldname)
            k_D .= 0.
            setfield!(spc, fieldname, k_D)
        end
    end
    return spc
end

function isolate_pmoas!(spc::AbstractParams, pmoas::Vector{String}, z::Int64)::Nothing
    spc = isolate_pmoas(spc, pmoasm, z)
    return nothing
end

function isolate_pmoas!(spc::AbstractParams, pmoas::Vector{String})::Nothing
    spc = isolate_pmoas(spc, pmoas)
    return nothing
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



@enum PMoA h G M A R # these are the valid PMoAs

"""
Apply zoom factor from reference species `deb_ref` to species of interest `spc`, 
given observed maximum structural mass `S_max` of the species of interest. \n

- `spc::AbstractParams`: default spc parameters of the species of interest
- `S_max::Float64`: observed maximum structural mass of the species of interest

kwargs:
- `deb_ref::AbstractParams`: spc parameters of the reference species, with default `SpeciesParams()`.
- `apply_to_covariates::Vector{Symbol}`: spc parameters which are assumed to covary with maximum structural mass. Zoom factor will be applied to these, assuming linear scaling.

"""
function zoom!(
    spc::AbstractParams,
    S_max::Float64;
    deb_ref::AbstractParams = SpeciesParams(),
    covariates::Vector{Symbol} = [:X_emb_int, :H_p, :K_X]
    )
    # TODO: update to exploit propagate_zoom field?
    S_max_ref = calc_S_max(deb_ref) # calculate maximum strucutral mass of the reference species 
    Z = S_max / S_max_ref # calculate the zoom factor 
    spc.Idot_max_rel *= Z^(1/3) # S_max scales with Idot_max_rel^3, so Z^(1/3) has to be paplied

    for covar in covariates # iterate over covariates
        previous_val = getproperty(spc, covar) # get the original value
        zoomed_val = previous_val * Z # apply the zoom factor 
        setproperty!(spc, covariate, zoomed_val) # update value in the parameter set
    end
end

"""
Set parameters with common prefix to the value of a reference parameter. 
E.g. 

```
set_equal!(spc, :Idot_max_rel, :lrv)
```

sets all parameters whose names start with `Idot_max_rel` equal to `Idot_max_rel_lrv`.
"""
function set_equal!(spc::AbstractParams, prefix::SS, ref_suffix::SS) where SS <: Union{Symbol,String}
    prefix = String(prefix)
    ref_suffix = String(ref_suffix)

    ref_paramname = prefix*"_"*ref_suffix
    ref_param = getproperty(spc, Symbol(ref_paramname))
    paramnames = String[String.(fieldnames(typeof(spc)))...]
    filter!(x -> occursin(String(prefix), x), paramnames)
    filter!(x -> x != String(ref_paramname), paramnames)

    for paramname in paramnames
        setproperty!(spc, Symbol(paramname), ref_param)
    end
end

# 
#import Base.setproperty!
"""
Throw warning if value assignemnt is done directly for DEB parameters. 
"""
#function setproperty!(spc::AbstractSpeciesParams, name::Symbol, value)
#    @warn("Direct assignment is discouraged for $(typeof(obj)), use assign!() instead")
#    setfield!(spc, name, value)
#end

