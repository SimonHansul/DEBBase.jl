"""
    agent_variability!(p::Ref{AbstractParams})
Induce agent variability in DEB parameters via zoom factor `Z`. 
`Z` is sampled from the corresponding distribution given in `p` and assumed to represent a ratio between maximum structurel *masses* (not lengths), 
so that the surface area-specific ingestion rate `Idot_max_rel` scales with `Z^(1/3)` and parameters which represent masses or energy pools scales with `Z`.
"""
function agent_variability!(pown::ComponentVector, pcmn::Ref{A}) where A <: AbstractParamCollection
    pown.Z = rand(pcmn.x.deb.Z) # sample zoom factor Z for agent a from distribution
    pown.Idot_max_rel = pcmn.x.deb.Idot_max_rel * pown.Z^(1/3) # Z is always applied to Idot_max_rel
    pown.Idot_max_rel_emb = pcmn.x.deb.Idot_max_rel_emb * pown.Z^(1/3) #, including the value for embryos

    for param in fieldnames(typeof(pcmn.x.deb.propagate_zoom)) # iterate over other parameters which may be affected by Z
        if getproperty(pcmn.x.deb.propagate_zoom, param) # check whether propagation of Z should occur for this parameter
            setproperty!(pcmn.x.deb, param, getproperty(pcmn.x.deb, param) * pown.Z) # assign the agent value by adjusting the hyperparameter
        else # if Z should not be propagated to this parameter, 
            setproperty!(pown, param, getproperty(pcmn.x.deb, param)) # set the agent-specific value equal to the population mean
        end
    end
end


"""
    initialize_statevars(pcmn::Ref{A}, pown::ComponentVector{Float64})::ComponentArray where A <: AbstractParamCollection
Initialize the component vector of state variables, `u`, based on common parameters `pcmn` and agent parameters `pown`.
"""
function initialize_statevars(pcmn::Ref{A}, pown::ComponentVector{Float64})::ComponentArray where A <: AbstractParamCollection
    return ComponentArray( # initial states
        X_p = Float64(pcmn.x.glb.Xdot_in), # initial resource abundance equal to influx rate
        X_emb = Float64(pown.X_emb_int), # initial mass of vitellus
        S = Float64(pown.X_emb_int * 0.001), # initial structure is a small fraction of initial reserve // mass of vitellus
        H = Float64(0.), # maturity
        H_b = 0., # maturity at birth (will be derived from model output)
        R = Float64(0.), # reproduction buffer
        I_emb = Float64(0.), # ingestion from vitellus
        I_p = Float64(0.), # ingestion from external food resource
        I = Float64(0.), # total ingestion
        A = Float64(0.), # assimilation
        M = Float64(0.), # somatic maintenance
        J = Float64(0.), # maturity maintenance 
        C_W = (pcmn.x.glb.C_W), # external stressor concentrations
        D_G = MVector{length(pcmn.x.deb.k_D_G), Float64}(zeros(length(pcmn.x.deb.k_D_G))), # scaled damage | growth efficiency
        D_M = MVector{length(pcmn.x.deb.k_D_M), Float64}(zeros(length(pcmn.x.deb.k_D_M))), # scaled damage | maintenance costs 
        D_A = MVector{length(pcmn.x.deb.k_D_A), Float64}(zeros(length(pcmn.x.deb.k_D_A))), # scaled damage | assimilation efficiency
        D_R = MVector{length(pcmn.x.deb.k_D_R), Float64}(zeros(length(pcmn.x.deb.k_D_R))), # scaled damage | reproduction efficiency
        D_h = MVector{length(pcmn.x.deb.k_D_h), Float64}(zeros(length(pcmn.x.deb.k_D_h))), # scaled damage | hazard rate
        y_G = Float64(1.), # relative response | growth efficiency
        y_M = Float64(1.), # relative response | maintenance costs 
        y_A = Float64(1.), # relative response | assimilation efficiency
        y_R = Float64(1.), # relative response | reproduction efficiency
        h_z = Float64(0.), # hazard rate | chemical stressors
        h_S = Float64(0.)  # hazard rate | starvation
    )
end

"""
    initialize_pown()::ComponentVector{Float64}
Initialize the agent-specific ("own") parameters. 
The parameters which can be agent-specific are predefined:
    - `Z`: zoom factor
    - `Idot_max_rel`
    - `Idot_max_rel_emb`
    - `X_emb_int`
    - `H_p`
    - `K_X`
"""
function initialize_pown()::ComponentVector{Float64}
    return ComponentVector{Float64}(
        Z = 1., 
        Idot_max_rel = 1e-310, 
        Idot_max_rel_emb = 1e-310, 
        X_emb_int = 1e-310, 
        H_p = 1e-310, 
        K_X = 1e-310
        )
end

"""
    simulator(
        pcmn::Ref{BaseParamCollection}; 
        saveat = 1,
        abstol = 1e-10, 
        reltol = 1e-10,
        kwargs...
        )::DataFrame

Run the DEBBase model from a reference to a `BaseParamCollection`  instance. <br>
These are the common parameters `pcmn`. The agent-specific parameters `pown` are initialized by the simulator. <br>
Additional kwargs are passed on to `DifferentialEquations.solve()`.
"""
function simulator(
    pcmn::Ref{BaseParamCollection}; 
    saveat = 1,
    abstol = 1e-10, 
    reltol = 1e-10,
    kwargs...
    )::DataFrame

    assert!(pcmn)
    pown = initialize_pown()
    agent_variability!(pown, pcmn)
    u = initialize_statevars(pcmn, pown)
    prob = ODEProblem(DEB!, u, (0, pcmn.x.glb.t_max), (pcmn, pown)) # define the problem
    sol = solve(prob, Tsit5(); saveat = saveat, abstol = abstol, reltol = reltol, kwargs...) # get solution to the IVP
    simout = sol_to_df(sol) # convert solution to dataframe
  
    return simout
end

"""
Run the DEBBase model from a `BaseParamCollection` instance. 
Creates a reference and runs the model from the reference. <br>
This is most convenient, but if many simulations are run (e.g. MC simulations, parameter sweeps), it less efficient than creating a reference, 
mutating the referenced object and running calling  `simulator` on the reference.
$(TYPEDSIGNATURES)
"""
function simulator(
    pcmn::BaseParamCollection;
    saveat = 1,
    abstol = 1e-10, 
    reltol = 1e-10,
    kwargs...
    )
    pref = Ref(pcmn) # create a reference to the parameter object
    return simulator(pref; saveat = saveat, abstol = abstol, reltol = reltol, kwargs...) # run simulation
end


"""
    @replicates(simcall::Expr, nreps::Int64) 

Perform replicated runs of `simcall`, where `simcall` is a call to a simulator function. 

Example:

    deb = DEBBaseParams(Z = Truncated(Normal(1, 0.1), 0, Inf)) # initialize default parameters with variable zoom factor
    yhat = @replicates DEBBase.simulator(BaseParamCollection(deb = deb))) 10 # execute replicated runs to simulator

In this case, `yhat` will contain the output of 10 replicated simulations. For each replicate, the zoom factor is sampled from a truncated Normal distribution. 
`yhat` contains an additional column `replicate`.
"""
macro replicates(simcall::Expr, nreps::Int64)
    quote
        yhat = DataFrame()

        for replicate in 1:$nreps
            yhat_i = $(esc(simcall))
            yhat_i[!,:replicate] .= replicate
            append!(yhat, yhat_i)
        end
        yhat
    end
end

