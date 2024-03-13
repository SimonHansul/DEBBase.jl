
@enum PMoA h G M A R # these are the valid PMoAs

"""
Apply zoom factor from reference species `deb_ref` to species of interest `deb`, 
given observed maximum structural mass `S_max` of the species of interest. \n

- `deb::AbstractParams`: default DEB parameters of the species of interest
- `S_max::Float64`: observed maximum structural mass of the species of interest

kwargs:
- `deb_ref::AbstractParams`: DEB parameters of the reference species, with default `DEBBaseParams()`.
- `apply_to_covariates::Vector{Symbol}`: DEB parameters which are assumed to covary with maximum structural mass. Zoom factor will be applied to these, assuming linear scaling.

"""
function zoom!(
    deb::AbstractParams,
    S_max::Float64;
    deb_ref::AbstractParams = DEBBaseParams(),
    covariates::Vector{Symbol} = [:X_emb_int, :H_p, :K_X]
    )
    S_max_ref = calc_S_max(deb_ref) # calculate maximum strucutral mass of the reference species 
    Z = S_max / S_max_ref # calculate the zoom factor 
    deb.Idot_max_rel *= Z^(1/3) # S_max scales with Idot_max_rel^3, so Z^(1/3) has to be paplied

    for covar in covariates # iterate over covariates
        previous_val = getproperty(deb, covar) # get the original value
        zoomed_val = previous_val * Z # apply the zoom factor 
        setproperty!(deb, covariate, zoomed_val) # update value in the parameter set
    end
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

    ref_paramname = prefix*"_"*ref_suffix
    ref_param = getproperty(deb, Symbol(ref_paramname))
    paramnames = String[String.(fieldnames(typeof(deb)))...]
    filter!(x -> occursin(String(prefix), x), paramnames)
    filter!(x -> x != String(ref_paramname), paramnames)

    for paramname in paramnames
        setproperty!(deb, Symbol(paramname), ref_param)
    end
end