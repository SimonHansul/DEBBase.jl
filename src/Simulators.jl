"""
Initialize the component vector of state variables, `u`, based on model parameters `p`.
$(TYPEDSIGNATURES)
"""
function initialize_statevars(p::AbstractParamCollection)
    return ComponentArray( # initial states
        X_p = p.glb.Xdot_in, # initial resource abundance equal to influx rate
        X_emb = p.deb.X_emb_int, # initial mass of vitellus
        S = p.deb.X_emb_int * 0.01, # initial structure is a small fraction of initial reserve // mass of vitellus
        H = 0., # maturity
        R = 0., # reproduction buffer
        life_stage = 0,
        I = 0.,
        A = 0.,
        M = 0.,
        J = 0.,
        D_G = zeros(length(p.deb.k_D_G)),
        D_M = zeros(length(p.deb.k_D_M)),
        D_A = zeros(length(p.deb.k_D_A)),
        D_R = zeros(length(p.deb.k_D_R)),
        D_h = zeros(length(p.deb.k_D_h)),
        C_W = p.glb.C_W,
        y_G = 1.,
        y_M = 1.,
        y_A = 1.,
        y_R = 1.,
        h_z = 0.,
        h_S = 0.
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
    
    # TODO: flatten D-Matrix into multiple columns in output dataframe
    sol = solve(prob, Tsit5(), reltol = 1e-6, abstol = 1e-10) # get solution to the IVP
    simout = sol_to_df(sol) # convert solution to dataframe
  
    return simout
end