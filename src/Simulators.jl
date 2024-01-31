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
        life_stage = Float64(1.), # life stage 
        I_emb = Float64(0.), # uptake from vitellus
        I_p = Float64(0.), # uptake from external food resource
        I = Float64(0.), # total uptake
        A = Float64(0.), # assimilation
        M = Float64(0.), # somatic maintenance
        J = Float64(0.), # maturity maintenance 
        C_W = (p.glb.C_W), # external stressor concentrations
        D_G = MVector{length(p.deb.k_D_G), Float64}(zeros(length(p.deb.k_D_G))), # scaled damage | growth efficiency
        D_M = MVector{length(p.deb.k_D_M), Float64}(zeros(length(p.deb.k_D_M))), # scaled damage | maintenance costs 
        D_A = MVector{length(p.deb.k_D_A), Float64}(zeros(length(p.deb.k_D_A))), # scaled damage | assimilation efficiency
        D_R = MVector{length(p.deb.k_D_R), Float64}(zeros(length(p.deb.k_D_R))), # scaled damage | reproduction efficiency
        D_h = MVector{length(p.deb.k_D_h), Float64}(zeros(length(p.deb.k_D_h))), # scaled damage | hazard rate
        y_G = Float64(1.), # relative response | growth efficiency
        y_M = Float64(1.), # relative response | maintenance costs 
        y_A = Float64(1.), # relative response | assimilation efficiency
        y_R = Float64(1.), # relative response | reproduction efficiency
        h_z = Float64(0.), # hazard rate | chemical stressors
        h_S = Float64(0.)  # hazard rate | starvation
    )
end

"""
Run any DEBBase-compatible model. 
$(TYPEDSIGNATURES)
"""
function simulator(
    p::BaseParamCollection; 
    saveat = 1,
    dt = 1/24
    )

    assert!(p)
    u = initialize_statevars(p)
    prob = ODEProblem(DEB!, u, (0, p.glb.t_max), p) # define the problem
    sol = solve(prob, Euler(), reltol = 1e-6, dt = dt, saveat = saveat) # get solution to the IVP
    #sol = solve(prob, Tsit5()) # get solution to the IVP
    simout = sol_to_df(sol) # convert solution to dataframe
  
    return simout
end