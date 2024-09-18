const X_EMB_INT_REL = 0.001 

function initialize_agent_statevars(p::Union{NamedTuple,AbstractParamCollection})
    ComponentVector(
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

function initalize_global_statevars(p::Union{NamedTuple,AbstractParamCollection})
    ComponentArray( # initial states
        X_p = p.glb.Xdot_in, # initial resource abundance equal to influx rate
        C_W = p.glb.C_W # external stressor concentrations
    )

end

"""
    initialize_statevars(p::AbstractParamCollection, pindt::ComponentVector{Float64})::ComponentArray

For initialization of ODE simulator, initialize the component vector of state variables, `u`, based on common oaraeter collection `p`.
"""
function initialize_statevars(p::Union{NamedTuple,AbstractParamCollection})::ComponentArray 
    vcat(
        initalize_global_statevars(p),
        initialize_agent_statevars(p)
    )
end