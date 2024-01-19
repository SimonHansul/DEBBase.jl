"""
Initialize the component vector of state variables, `u`, based on model parameters `p`.
$(TYPEDSIGNATURES)
"""
function initialize_statevars(p::BaseParamCollection)::ComponentArray
    return ComponentArray( # initial states
        X_p = Float64(p.glb.Xdot_in), # initial resource abundance equal to influx rate
        X_emb = Float64(p.deb.X_emb_int), # initial mass of vitellus
        S = Float64(p.deb.X_emb_int * 0.01), # initial structure is a small fraction of initial reserve // mass of vitellus
        H = Float64(0.), # maturity
        R = Float64(0.), # reproduction buffer
        life_stage = 1., # life stage 
        I_emb = 0., # uptake from vitellus
        I_p = 0., # uptake from external food resource
        I = 0., # total uptake
        A = 0., # assimilation
        M = 0., # somatic maintenance
        J = 0., # maturity maintenance 
        C_W = p.glb.C_W, # external stressor concentrations
        D_G = zeros(length(p.deb.k_D_G)), # scaled damage | growth efficiency
        D_M = zeros(length(p.deb.k_D_M)), # scaled damage | maintenance costs 
        D_A = zeros(length(p.deb.k_D_A)), # scaled damage | assimilation efficiency
        D_R = zeros(length(p.deb.k_D_R)), # scaled damage | reproduction efficiency
        D_h = zeros(length(p.deb.k_D_h)), # scaled damage | hazard rate
        y_G = 1., # relative response | growth efficiency
        y_M = 1., # relative response | maintenance costs 
        y_A = 1., # relative response | assimilation efficiency
        y_R = 1., # relative response | reproduction efficiency
        h_z = 0., # hazard rate | chemical stressors
        h_S = 0.  # hazard rate | starvation
    )
end

"""
Run any DEBBase-compatible model. 
$(TYPEDSIGNATURES)
"""
function simulator(
    p::BaseParamCollection
    )

    assert!(p)
    u = initialize_statevars(p)
    prob = ODEProblem(DEB!, u, (0, p.glb.t_max), p) # define the problem
    sol = solve(prob, Euler(), reltol = 1e-6, dt = 1/24) # get solution to the IVP
    simout = sol_to_df(sol) # convert solution to dataframe
  
    return simout
end