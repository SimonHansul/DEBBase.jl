
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