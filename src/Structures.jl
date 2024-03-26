@with_kw mutable struct GlobalBaseParams <: AbstractParams
    t_max::Float64 = 21.
    Xdot_in::Float64 = 1200. # set to be a little above the absolute maximum ingestion rate according to default DEBBaseParams
    V_patch::Float64 = 0.05
    C_W::Vector{Float64} = [0.]
    units::NamedTuple = (time = "d", mass = "mug C", volume = "L")
end

"""
DEBBase Parameters. <br>
Default values are for Daphnia magna and Azoxystrobin.
"""
@with_kw mutable struct DEBBaseParams <: AbstractParams
    Z::Distribution = Truncated(Normal(1., 1.), 0, Inf) # Individual variability is accounted for in the zoom factor. This can be set to a Dirac distribution if a zoom factor should be applied without introducing individual variability.
    propagate_zoom::Vector{Symbol} = [:X_emb_int, :H_p, :K_X] # Parameters to which Z will be propagated. Z is *always* applied to `Idot_max_rel` (with appropriate scaling).
    X_emb_int_mean::Float64 = 19.42 # Population mean of the initial vitellus mass
    X_emb_int::Float64 = NaN # Individual-specific initial vitellus (≈ yolk) mass
    K_X_mean::Float64 = 1. # Population mean of half-saturation constant for food uptake
    K_X::Float64 = 1. # Individual half-saturation constant for food uptake
    Idot_max_rel_mean::Float64 = 22.9 # Population mean if the maximum size-specific ingestion rate
    Idot_max_rel::Float64 = NaN # Individual maximum size-specific ingestion rate
    Idot_max_rel_emb_mean::Float64 = 22.9 # Population mean of the size-specific embryonic ingestion rate
    Idot_max_rel_emb::Float64 = 22.9 # Individual size-specific embryonic ingestion rate
    kappa::Float64 = 0.539 # Somatic allocation fraction
    eta_IA::Float64 = 0.33 # Assimilation efficiency
    eta_AS::Float64 = 0.8 # Growth efficiency
    eta_SA::Float64 = 0.8 # Shrinking efficiency
    eta_AR::Float64 = 0.95 # Reproduction efficiency
    k_M::Float64 = 0.59 # Somatic maintenance rate constant
    k_J::Float64 = 0.504 # Mautirty maintenance rate constant
    H_p_mean::Float64 = 100. # Population mean of maturity at puberty
    H_p::Float64 = NaN # Individual maturity at puberty
    
    k_D_G::Vector{Float64} = [0.] # Dominant rate constants | PMoA growth efficiency
    k_D_M::Vector{Float64} = [0.] # Dominant rate constants | PMoA maintenance costs
    k_D_A::Vector{Float64} = [0.] # Dominant rate constants | PMoA assimilation efficiency
    k_D_R::Vector{Float64} = [0.38] # Dominant rate constants | PMoA reproduction efficiency
    k_D_h::Vector{Float64} = [0.] # Dominant rate constants | PMoA hazard rate
    
    drc_functs_G::Vector{Function} = [LL2] # Dose-response functions | PMoA growth efficiency
    drc_functs_M::Vector{Function} = [LL2M] # Dose-response functions | PMoA maintenance costs
    drc_functs_A::Vector{Function} = [LL2] # Dose-response functions | PMoA assimilation efficiency
    drc_functs_R::Vector{Function} = [LL2] # Dose-response functions | PMoA reproduction efficiency
    drc_functs_h::Vector{Function} = [LL2h] # Dose-response functions | PMoA hazard rate

    drc_params_G::Vector{NTuple} = [(1e10, 1e10)] # Dose-response parameters | PMoA growth efficiency
    drc_params_M::Vector{NTuple} = [(1e10, 1e10)] # Dose-response parameters | PMoA maintenance costs
    drc_params_A::Vector{NTuple} = [(1e10, 1e10)] # Dose-response parameters | assimilation efficiency
    drc_params_R::Vector{NTuple} = [(167., 0.93)] # Dose-response parameters | PMoA reproduction efficiency
    drc_params_h::Vector{NTuple} = [(1e10, 1e10)] # Dose-response parameters | PMoA hazard rate
end

@with_kw mutable struct BaseParamCollection <: AbstractParamCollection
    glb::AbstractParams = GlobalBaseParams()
    deb::AbstractParams = DEBBaseParams()
end



