abstract type Params end

@with_kw mutable struct GlobalParams <: Params
    t_max::Float64 = 21.
    Xdot_in::Float64 = 10.
end

@with_kw mutable struct DEBParams <: Params
    X_emb_int::Float64 = 0.001
    K_X::Float64 = 1.
    Idot_max_rel::Float64 = 0.75
    kappa::Float64 = 0.75
    eta_IA::Float64 = 0.8
    eta_AS::Float64 = 0.8
    eta_SA::Float64 = 0.8
    eta_AR::Float64 = 0.95
    k_M::Float64 = 0.1
    k_J::Float64 = 0.08
    H_p::Float64 = 2.
end