# Abstract types are defined so that dispatch works as expected across model versions and Parameter / Statevar types 
# For example, we can define a Vector{AbstractParams} to collect Global and DEB Params. 

abstract type AbstractOrganism end
abstract type AbstractParams end
abstract type AbstractStatevars end

@with_kw mutable struct GlobalBaseParams <: AbstractParams
    t_max::Float64 = 21.
    Xdot_in::Float64 = 1200. # set to be a little above the absolute maximum ingestion rate according to default DEBBaseParams
    V_patch = 0.05
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
end

@with_kw mutable struct DEBBaseStatevars <: AbstractStatevars
    X_emb::Float64
    S::Float64
    H::Float64 = 0.
    R::Float64 = 0.
    life_stage::Float64 = 1.
    Idot::Float64 = 0.
end

@with_kw mutable struct DEBBaseOrganism <: AbstractOrganism
    debparams::DEBBaseParams = DEBBaseParams()
    statevars::DEBBaseStatevars = DEBBaseStatevars(
        X_emb = DEBBaseParams().X_emb_int,
        S = DEBBaseParams().X_emb_int * 0.01
        )
end