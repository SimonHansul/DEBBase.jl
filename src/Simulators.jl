"""
Initialize the component vector of state variables, `u`, based on model parameters `p`.
$(TYPEDSIGNATURES)
"""
function initialize_statevars(p::Ref{A})::ComponentArray where A <: AbstractParamCollection
    return ComponentArray( # initial states
        X_p = Float64(p.x.glb.Xdot_in), # initial resource abundance equal to influx rate
        X_emb = Float64(p.x.deb.X_emb_int), # initial mass of vitellus
        S = Float64(p.x.deb.X_emb_int * 0.001), # initial structure is a small fraction of initial reserve // mass of vitellus
        H = Float64(0.), # maturity
        H_b = 0., # maturity at birth (will be derived from model output)
        R = Float64(0.), # reproduction buffer
        I_emb = Float64(0.), # uptake from vitellu
        I_p = Float64(0.), # uptake from external food resource
        I = Float64(0.), # total uptake
        A = Float64(0.), # assimilation
        M = Float64(0.), # somatic maintenance
        J = Float64(0.), # maturity maintenance 
        C_W = (p.x.glb.C_W), # external stressor concentrations
        D_G = MVector{length(p.x.deb.k_D_G), Float64}(zeros(length(p.x.deb.k_D_G))), # scaled damage | growth efficiency
        D_M = MVector{length(p.x.deb.k_D_M), Float64}(zeros(length(p.x.deb.k_D_M))), # scaled damage | maintenance costs 
        D_A = MVector{length(p.x.deb.k_D_A), Float64}(zeros(length(p.x.deb.k_D_A))), # scaled damage | assimilation efficiency
        D_R = MVector{length(p.x.deb.k_D_R), Float64}(zeros(length(p.x.deb.k_D_R))), # scaled damage | reproduction efficiency
        D_h = MVector{length(p.x.deb.k_D_h), Float64}(zeros(length(p.x.deb.k_D_h))), # scaled damage | hazard rate
        y_G = Float64(1.), # relative response | growth efficiency
        y_M = Float64(1.), # relative response | maintenance costs 
        y_A = Float64(1.), # relative response | assimilation efficiency
        y_R = Float64(1.), # relative response | reproduction efficiency
        h_z = Float64(0.), # hazard rate | chemical stressors
        h_S = Float64(0.)  # hazard rate | starvation
    )
end


"""
Run the DEBBase model from a reference to a `BaseParamCollection`  instance.
$(TYPEDSIGNATURES)
"""
function simulator(
    p::Ref{BaseParamCollection}; 
    alg,
    saveat = 1,
    abstol = 1e-10, 
    reltol = 1e-10,
    kwargs...
    )

    assert!(p)
    u = initialize_statevars(p)
    prob = ODEProblem(DEB!, u, (0, p.glb.t_max), p) # define the problem
    sol = solve(prob, alg, reltol = reltol, saveat = saveat; kwargs...) # get solution to the IVP
    simout = sol_to_df(sol) # convert solution to dataframe
  
    return simout
end

"""
Run the DEBBase model from a `BaseParamCollection` instance. 
Creates a reference and runs the model from the reference. <br>
This is most convenient, but if many simulations are run (e.g. MC simulations, parameter sweeps), 
this is less efficient than creating a reference, 
mutating the referenced object and running calling  `simulator` on the reference.
$(TYPEDSIGNATURES)
"""
function simulator(
    p::BaseParamCollection;
    saveat = 1,
    abstol = 1e-10, 
    reltol = 1e-10,
    kwargs...
    )
    pref = Ref(p) # create a reference to the parameter object
    return simulator(pref; saveat = saveat, abstol = abstol, reltol = reltol, kwargs...) # run simulation
end

