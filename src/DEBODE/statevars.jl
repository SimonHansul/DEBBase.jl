const X_EMB_INT_REL = 0.001 


"""
    initialize_statevars(theta::Union{AbstractParamCollection,NamedTuple})::ComponentArray 

For initialization of ODE simulator, initialize the component vector of state variables, `u`, based on common parameter collection `theta`.
"""
function initialize_statevars(p::Union{AbstractParamCollection,NamedTuple})::ComponentArray 
    return ComponentArray( # initial states
        X_p = p.glb.Xdot_in, # initial resource abundance equal to influx rate
        C_W = p.glb.C_W, # external stressor concentrations
        T = p.glb.T, # ambient temperature

        embryo = 1, # life stage indicators - determined in callbacks
        juvenile = 0,
        adult = 0,

        X_emb = p.agn.X_emb_int_0, # initial mass of vitellus
        S = p.agn.X_emb_int_0 * X_EMB_INT_REL, # initial structure is a small fraction of initial reserve // mass of vitellus
        S_max_hist = p.agn.X_emb_int_0 * X_EMB_INT_REL, # initial reference structure
        H = 0, # maturity
        H_b = 0, # maturity at birth (will be derived from model output)
        R = 0, # reproduction buffer
        f_X = 1, # scaled functional response 
        I_emb = 0, # ingestion from vitellus
        I_p = 0, # ingestion from external food resource
        I = 0, # total ingestion
        A = 0, # assimilation
        M = 0, # somatic maintenance
        J = 0, # maturity maintenance 
        Q = 0, # cumulative dissipation flux

        X_emb_int = p.agn.X_emb_int_0,
        Idot_max_rel = p.agn.Idot_max_rel_0, 
        Idot_max_rel_emb = p.agn.Idot_max_rel_emb_0, 
        K_X = p.agn.K_X_0,
        eta_AS = p.spc.eta_AS_0, # current growth efficiency
        kappa = p.spc.kappa_0, 
        k_M = p.spc.k_M_0, # current maintenance rate constant
        k_J = p.spc.k_J_0, # current maturity rate constant
        eta_IA = p.spc.eta_IA_0, # current assimilation efficiency
        eta_AR = p.spc.eta_AR_0, # current reproduction efficiency
        # H_p == H_p_0 in the base model, but we need to access it in lifestage transition callbacks
        # #TODO: how to access a parameter within a nested struct in a callback?
        H_p = p.agn.H_p_0, # current maturity threshold at puberty; 

        D_G = MVector{length(p.spc.k_D_G), Float64}(zeros(length(p.spc.k_D_G))), # scaled damage | growth efficiency
        D_M = MVector{length(p.spc.k_D_M), Float64}(zeros(length(p.spc.k_D_M))), # scaled damage | maintenance costs 
        D_A = MVector{length(p.spc.k_D_A), Float64}(zeros(length(p.spc.k_D_A))), # scaled damage | assimilation efficiency
        D_R = MVector{length(p.spc.k_D_R), Float64}(zeros(length(p.spc.k_D_R))), # scaled damage | reproduction efficiency
        D_h = MVector{length(p.spc.k_D_h), Float64}(zeros(length(p.spc.k_D_h))), # scaled damage | hazard rate

        y_G = 1., # relative response | growth efficiency
        y_M = 1., # relative response | maintenance costs 
        y_A = 1., # relative response | assimilation efficiency
        y_R = 1., # relative response | reproduction efficiency
        y_T = 1., # relative response | temperature effects
        h_z = 0. # hazard rate | chemical stressors
    )
end