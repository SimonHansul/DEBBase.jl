
"""
Definition of the DEB model for use with ODE solver (du = derivatives, u = state variables, p = parameters, t = time), 
including stressor effects.
$(TYPEDSIGNATURES)
"""
function DEB!(du, u, p, t)
    let y_G, y_M, y_A, y_R, y_E, y_K # define local variables
        Xdot_in, # unpack parameters
        volume,
        C_1,
        C_2,
        X_emb_int,
        K_X,
        Idot_max_rel_0_lrv,
        Idot_max_rel_0_ad,
        kappa_0_lrv,
        kappa_0_ad,
        eta_IA_0_lrv,
        eta_IA_0_ad,
        eta_AS_0_lrv,
        eta_AS_0_ad,
        eta_SA,
        eta_AR_0,
        k_M_0_lrv,
        k_M_0_mt,
        k_M_0_ad,
        k_J_0,
        k_C,
        H_42,
        H_46,
        H_p,
        k_d_h_1,
        k_d_G_1,
        k_d_M_1,
        k_d_A_1,
        k_d_R_1,
        k_d_E_1,
        k_d_K_1,
        k_d_h_2,
        k_d_G_2,
        k_d_M_2,
        k_d_A_2,
        k_d_R_2,
        k_d_E_2,
        k_d_K_2,
        drcmodel_h_1,
        drcmodel_G_1,
        drcmodel_M_1,
        drcmodel_A_1,
        drcmodel_R_1,
        drcmodel_E_1,
        drcmodel_K_1,
        drcmodel_h_2,
        drcmodel_G_2,
        drcmodel_M_2,
        drcmodel_A_2,
        drcmodel_R_2,
        drcmodel_E_2,
        drcmodel_K_2,
        drcparams_h_1,
        drcparams_G_1,
        drcparams_M_1,
        drcparams_A_1,
        drcparams_R_1,
        drcparams_E_1,
        drcparams_K_1,
        drcparams_h_2,
        drcparams_G_2,
        drcparams_M_2,
        drcparams_A_2,
        drcparams_R_2,
        drcparams_E_2,
        drcparams_K_2,
        mixture_model,
        metamorphic_transition_function = p
        
        # unpack state variables
        X_p,
        X_emb,
        S,
        H,
        R,
        D_h_1,
        D_G_1,
        D_M_1,
        D_A_1,
        D_R_1,
        D_E_1,
        D_K_1,
        D_h_2,
        D_G_2,
        D_M_2,
        D_A_2,
        D_R_2,
        D_E_2,
        D_K_2,
        S_z = u

        D_h_1 = max(0, D_h_1)
        D_G_1 = max(0, D_G_1)
        D_M_1 = max(0, D_M_1)
        D_A_1 = max(0, D_A_1)
        D_R_1 = max(0, D_R_1)
        D_E_1 = max(0, D_E_1)
        D_K_1 = max(0, D_K_1)

        D_h_2 = max(0, D_h_2)
        D_G_2 = max(0, D_G_2)
        D_M_2 = max(0, D_M_2)
        D_A_2 = max(0, D_A_2)
        D_R_2 = max(0, D_R_2)
        D_E_2 = max(0, D_E_2)
        D_K_2 = max(0, D_K_2)

        S = max(0, S) # control for negative values

        life_stage = determine_life_stage(H, H_42, H_46, H_p, X_emb, X_emb_int)
        Idot_max_rel_0, eta_IA_0, k_M_0, kappa_0, eta_AS_0 = life_stage_effects(
            H_42,
            H_46,
            H,
            metamorphic_transition_function,
            juvenile(life_stage),
            adult(life_stage),
            metamorph(life_stage),
            Idot_max_rel_0_lrv,
            Idot_max_rel_0_ad,
            eta_IA_0_lrv,
            eta_IA_0_ad,
            k_M_0_lrv,
            k_M_0_mt,
            k_M_0_ad,
            kappa_0_lrv,
            kappa_0_ad,
            eta_AS_0_lrv,
            eta_AS_0_ad,
        )

        # calculate intermediate variables for TKTD submodel
        h_1 = drcmodel_h_1(D_h_1, drcparams_h_1) # calculate individual hazard rates
        h_2 = drcmodel_h_2(D_h_2, drcparams_h_2)

        dD_h_1 = minimaltk(D_h_1, C_1, k_d_h_1)
        dD_h_2 = minimaltk(D_h_2, C_2, k_d_h_2)
        dD_G_1 = minimaltk(D_G_1, C_1, k_d_G_1)
        dD_G_2 = minimaltk(D_G_2, C_2, k_d_G_2)
        dD_M_1 = minimaltk(D_M_1, C_1, k_d_M_1)
        dD_M_2 = minimaltk(D_M_2, C_2, k_d_M_2)
        dD_A_1 = minimaltk(D_A_1, C_1, k_d_A_1)
        dD_A_2 = minimaltk(D_A_2, C_2, k_d_A_2)
        dD_R_1 = minimaltk(D_R_1, C_1, k_d_R_1)
        dD_R_2 = minimaltk(D_R_2, C_2, k_d_R_2)
        dD_E_1 = minimaltk(D_E_1, C_1, k_d_E_1)
        dD_E_2 = minimaltk(D_E_2, C_2, k_d_E_2)
        dD_K_1 = minimaltk(D_K_1, C_1, k_d_K_1)
        dD_K_2 = minimaltk(D_K_2, C_2, k_d_K_2)

        y_G_1 = drcmodel_G_1(D_G_1, drcparams_G_1)
        y_G_2 = drcmodel_G_2(D_G_2, drcparams_G_2)
        y_M_1 = drcmodel_M_1(D_M_1, drcparams_M_1)
        y_M_2 = drcmodel_M_2(D_M_2, drcparams_M_2)

        y_A_1 = drcmodel_A_1(D_A_1, drcparams_A_1)
        y_A_2 = drcmodel_A_2(D_A_2, drcparams_A_2)
        y_R_1 = drcmodel_R_1(D_R_1, drcparams_R_1)
        y_R_2 = drcmodel_R_2(D_R_2, drcparams_R_2)
        y_E_1 = drcmodel_E_1(D_E_1, drcparams_E_1)
        y_E_2 = drcmodel_E_2(D_E_2, drcparams_E_2)
        y_K_1 = drcmodel_K_1(D_K_1, drcparams_K_1)
        y_K_2 = drcmodel_K_2(D_K_2, drcparams_K_2)

        if mixture_model == "independent_action"
            h = h_1 + h_2 # hazard rates are additive

            y_G = y_G_1 * y_G_2 # relative responses are multiplicative
            y_M = y_M_1 * y_M_2
            y_A = y_A_1 * y_A_2
            y_R = y_R_1 * y_R_2
            y_E = y_E_1 * y_E_2
            y_K = y_K_1 * y_K_2

        elseif mixture_model == "damage_addition"
            h = drcmodel_h_1(D_h_1 + D_h_2, drcparams_h_1)

            y_G = drcmodel_G_1(D_G_1 + D_G_2, drcparams_G_1)
            y_M = drcmodel_M_1(D_M_1 + D_M_2, drcparams_M_1)
            y_A = drcmodel_A_1(D_A_1 + D_A_2, drcparams_A_1)
            y_R = drcmodel_R_1(D_R_1 + D_R_2, drcparams_R_1)
            y_E = drcmodel_E_1(D_E_1 + D_E_2, drcparams_E_1)
            y_K = drcmodel_K_1(D_K_1 + D_K_2, drcparams_K_1)
        end

        # calculate survival probability assuming SD mechanism
        s_zdot = -h * S_z

        # apply stress to DEB processes
        eta_AS = eta_AS_0 * y_G
        k_M = k_M_0 * y_M
        k_J = k_J_0 * y_M
        eta_IA = eta_IA_0 * y_A
        eta_AR = eta_AR_0 * y_R
        # effects on X_emb int not relevant for pODE implementation
        kappa = kappa_0 # kappa, ingestion rate and searching rate not (yet) considered as PMoA
        Idot_max_rel = Idot_max_rel_0 
        
        idot = Idot(X_p, volume, X_emb, life_stage, S, Idot_max_rel, K_X)
        xdot = embryo(life_stage) ? Xdot_in : Xdot_in - idot
        xembdot = embryo(life_stage) ? -idot : 0.0
        cdot = Cdot(S, k_C, metamorph(life_stage))
        adot = idot * eta_IA + cdot * eta_SA
        mdot = Mdot(S, k_M)
        jdot = Jdot(H, k_J)
        sdot = Sdot(adot, mdot, kappa, eta_AS, eta_SA, cdot) 
        hdot = adult(life_stage) ? 0.0 : (1 - kappa) * adot - jdot
        rdot = adult(life_stage) ? eta_AR * (1 - kappa) * adot - jdot : 0.0

        # update du/dt
        du[01] = xdot
        du[02] = xembdot
        du[03] = sdot
        du[04] = hdot
        du[05] = rdot
        du[06] = dD_h_1
        du[07] = dD_G_1
        du[08] = dD_M_1
        du[09] = dD_A_1
        du[10] = dD_R_1
        du[11] = dD_E_1
        du[12] = dD_K_1
        du[13] = dD_h_2
        du[14] = dD_G_2
        du[15] = dD_M_2
        du[16] = dD_A_2
        du[17] = dD_R_2
        du[18] = dD_E_2
        du[19] = dD_K_2
        du[20] = s_zdot
    end # let
end
