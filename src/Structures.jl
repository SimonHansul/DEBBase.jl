# Abstract types are defined so that dispatch works as expected across model versions and Parameter / Statevar types 
# For example, we can define a Vector{AbstractParams} to collect Global and DEB Params. 

abstract type AbstractOrganism end
abstract type AbstractParams end
abstract type AbstractStatevars end

@with_kw mutable struct GlobalBaseParams <: AbstractParams
    t_max::Float64 = 21.
    Xdot_in::Float64 = 1200. # set to be a little above the absolute maximum ingestion rate according to default DEBBaseParams
    V_patch::Float64 = 0.05
    C_W::Vector{Float64} = [0., 0.]
    units::NamedTuple = (time = "d", mass = "mug", volume = "L")
end

@with_kw mutable struct GlobalBaseStatevars <: AbstractStatevars
    X_p::Float64 = 0.
end

"""
DEBBase Parameters with default values for Daphnia magna.
$(TYPEDSIGNATURES)
"""
@with_kw mutable struct DEBBaseParams <: AbstractParams
    # physiological baseline (DEB) parameters
    X_emb_int::Float64 = 19.42
    K_X::Float64 = 1.
    Idot_max_rel::Float64 = 22.9
    Idot_max_rel_emb::Float64 = 22.9
    kappa::Float64 = 0.539
    eta_IA::Float64 = 0.33
    eta_AS::Float64 = 0.8
    eta_SA::Float64 = 0.8
    eta_AR::Float64 = 0.95
    k_M::Float64 = 0.59
    k_J::Float64 = 0.
    H_p::Float64 = 100.
    
    # dominant rate constants
    k_D_G::Vector{Float64} = [1., 1.]
    k_D_M::Vector{Float64} = [1., 1.]
    k_D_A::Vector{Float64} = [1., 1.]
    k_D_R::Vector{Float64} = [1., 1.]
    k_D_h::Vector{Float64} = [1., 1.]
    
    # DRC functions per PMoA
    drc_functs_G::Vector{Function} = [LL2, LL2]
    drc_functs_M::Vector{Function} = [LL2M, LL2M]
    drc_functs_A::Vector{Function} = [LL2, LL2]
    drc_functs_R::Vector{Function} = [LL2, LL2]
    drc_functs_h::Vector{Function} = [LL2h, LL2h]

    # DRC parameters per PMoA
    drc_params_G::Vector{NTuple} = [(1., 2.), (1., 2.)]
    drc_params_M::Vector{NTuple} = [(1., 2.), (1., 2.)]
    drc_params_A::Vector{NTuple} = [(1., 2.), (1., 2.)]
    drc_params_R::Vector{NTuple} = [(1., 2.), (1., 2.)]
    drc_params_h::Vector{NTuple} = [(1., 2.), (1., 2.)]
end