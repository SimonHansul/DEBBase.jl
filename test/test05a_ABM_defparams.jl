# script to test the DEB-ABM on the default parameter set
using Pkg; Pkg.activate("test")
using Parameters
using NamedTupleTools

using Revise 
@time using DEBBase.DEBODE
using DEBBase.ParamStructs
using DEBBase.DoseResponse

@with_kw mutable struct GlobalABMParams <: ParamStructs.AbstractGlobalParams
    N0::Int64 = 1 #  initial number of individuals [#]
    t_max::Float64 = 21. # maximum simulation time [d]
    Xdot_in::Float64 = 1200. # nutrient influx set to be a little above the absolute maximum ingestion rate according to default SpeciesParams [μg C d-1]
    k_V::Float64 = 0. # chemostat dilution rate [d-1]
    V_patch::Float64 = 0.5 # volume of a patch (L) (or the entire similated environment) [L]
    C_W::Vector{Float64} = [0.] # external chemical concentrations [μg L-1]
end

@with_kw mutable struct SpeciesABMParams <: AbstractSpeciesParams
    Z::Distribution = Dirac(1.) # agent variability is accounted for in the zoom factor. This can be set to a Dirac distribution if a zoom factor should be applied without introducing agent variability.
    propagate_zoom::@NamedTuple{X_emb_int::Bool, H_p::Bool, K_X::Bool} = (X_emb_int = true, H_p = true, K_X = true) # Parameters to which Z will be propagated. Z is *always* applied to `Idot_max_rel` (with appropriate scaling).
    X_emb_int::Float64 = 19.42 # initial vitellus [μgC]
    K_X::Float64 = 1. # half-saturation constant for food uptake [μgC L-1]
    Idot_max_rel::Float64 = 22.9 # maximum size-specific ingestion rate [μgC μgC^-(2/3) d-1]
    Idot_max_rel_emb::Float64 = 22.9 # size-specific embryonic ingestion rate [μgC μgC^-(2/3) d-1]
    kappa::Float64 = 0.539 # somatic allocation fraction [-]
    eta_IA::Float64 = 0.33 # assimilation efficiency [-]
    eta_AS::Float64 = 0.8 # growth efficiency [-]
    eta_SA::Float64 = 0.8 # shrinking efficiency [-]
    eta_AR::Float64 = 0.95 # reproduction efficiency [-]
    k_M::Float64 = 0.59 # somatic maintenance rate constant [d^-1]
    k_J::Float64 = 0.504 # maturity maintenance rate constant [d^-1]
    H_p::Float64 = 100. # maturity at puberty [μgC]
    
    k_D_G::Vector{Float64} = [0.00] # toxicokinetic rate constants | PMoA growth efficiency
    k_D_M::Vector{Float64} = [0.00] # toxicokinetic rate constants | PMoA maintenance costs
    k_D_A::Vector{Float64} = [0.00] # toxicokinetic rate constants | PMoA assimilation efficiency
    k_D_R::Vector{Float64} = [0.38] # toxicokinetic rate constants | PMoA reproduction efficiency
    k_D_h::Vector{Float64} = [0.00] # toxicokinetic rate constants | PMoA hazard rate
    
    drc_functs_G::Vector{Function} = [LL2]  # dose-response functions | PMoA growth efficiency
    drc_functs_M::Vector{Function} = [LL2M] # dose-response functions | PMoA maintenance costs
    drc_functs_A::Vector{Function} = [LL2]  # dose-response functions | PMoA assimilation efficiency
    drc_functs_R::Vector{Function} = [LL2]  # dose-response functions | PMoA reproduction efficiency
    drc_functs_h::Vector{Function} = [LL2h] # dose-response functions | PMoA hazard rate

    e_G::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA growth efficiency
    e_M::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA maintenance costs
    e_A::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA assimilation efficiency
    e_R::Union{Nothing,Vector{Float64}} = [167.] # sensitivity parameters | PMoA reproduction efficiency
    e_h::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA hazard rate

    b_G::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA growth efficiency
    b_M::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA maintenance costs
    b_A::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA assimilation efficiency
    b_R::Union{Nothing,Vector{Float64}} = [0.93] # slope parameters | PMoA reproduction efficiency
    b_h::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA reproduction efficiency 
end

# first we extend the species parameters to include additional params needed for the ABM

glb = GlobalABMParams() |> ntfromstruct
spc_ode = SpeciesABMParams() |> ntfromstruct
agn = DEBODE.ODEAgentParams(DEBParamCollection())

spc = (; 
    spc_ode...,
    (
        a_max = 60., # maximum life span
        tau_R = 2. # reproduction interval
    )...
)






@with_kw mutable struct DEBAgent
    du::ComponentVector
    u::ComponentVector
    p::Union{AbstractParamCollection,NamedTuple}
    
end