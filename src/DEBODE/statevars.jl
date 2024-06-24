const X_EMB_INT_REL = 0.001 

"""
    initialize_statevars(p::AbstractParamCollection, pindt::ComponentVector{Float64})::ComponentArray

For initialization of ODE simulator, initialize the component vector of state variables, `u`, based on common oaraeter collection `theta`.
"""
function initialize_statevars(theta::AbstractParamCollection)::ComponentArray 
    return ComponentArray( # initial states
        X_p = theta.glb.Xdot_in, # initial resource abundance equal to influx rate
        C_W = theta.glb.C_W, # external stressor concentrations

        X_emb = theta.agn.X_emb_int, # initial mass of vitellus
        S = theta.agn.X_emb_int * X_EMB_INT_REL, # initial structure is a small fraction of initial reserve // mass of vitellus
        S_max_hist = theta.agn.X_emb_int * X_EMB_INT_REL, # initial reference structure
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

        D_G = MVector{length(theta.spc.k_D_G), Float64}(zeros(length(theta.spc.k_D_G))), # scaled damage | growth efficiency
        D_M = MVector{length(theta.spc.k_D_M), Float64}(zeros(length(theta.spc.k_D_M))), # scaled damage | maintenance costs 
        D_A = MVector{length(theta.spc.k_D_A), Float64}(zeros(length(theta.spc.k_D_A))), # scaled damage | assimilation efficiency
        D_R = MVector{length(theta.spc.k_D_R), Float64}(zeros(length(theta.spc.k_D_R))), # scaled damage | reproduction efficiency
        D_h = MVector{length(theta.spc.k_D_h), Float64}(zeros(length(theta.spc.k_D_h))), # scaled damage | hazard rate

        y_G = 1., # relative response | growth efficiency
        y_M = 1., # relative response | maintenance costs 
        y_A = 1., # relative response | assimilation efficiency
        y_R = 1., # relative response | reproduction efficiency
        h_z = 0. # hazard rate | chemical stressors
    )
end