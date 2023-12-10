@with_kw mutable struct DEBParams
    # TODO: 
    # Add xenopus params from IR2 as defaults
    # but adjust k_M_0_mt and maybe H_46
    X_emb_int::Float64
    K_X::Float64
    Idot_max_rel_0_lrv::Float64
    Idot_max_rel_0_ad::Float64 
    kappa_0_lrv::Float64
    kappa_0_ad::Float64
    eta_IA_0_lrv::Float64
    eta_IA_0_ad::Float64
    eta_AS_0_lrv::Float64
    eta_AS_0_ad::Float64
    eta_SA::Float64
    eta_AR_0::Float64
    k_M_0_lrv::Float64
    k_M_0_mt::Float64
    k_M_0_ad::Float64
    k_J_0::Float64
    k_C::Float64
    H_42::Float64
    H_46::Float64
    H_p::Float64
end