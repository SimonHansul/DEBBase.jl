


"""
    initialize_statevars(p::AbstractParamCollection, pagnt::ComponentVector{Float64})::ComponentArray 
Initialize the component vector of state variables, `u`, based on common parameters `p` and agent parameters `pagnt`.
"""
function initialize_statevars(p::AbstractParamCollection)::ComponentArray 
    return ComponentArray( # initial states
        X_p = Float64(p.glb.Xdot_in), # initial resource abundance equal to influx rate
        X_emb = Float64(p.agn.X_emb_int), # initial mass of vitellus
        S = Float64(p.agn.X_emb_int * 0.001), # initial structure is a small fraction of initial reserve // mass of vitellus
        H = Float64(0.), # maturity
        H_b = 0., # maturity at birth (will be derived from model output)
        R = Float64(0.), # reproduction buffer
        I_emb = Float64(0.), # ingestion from vitellus
        I_p = Float64(0.), # ingestion from external food resource
        I = Float64(0.), # total ingestion
        A = Float64(0.), # assimilation
        M = Float64(0.), # somatic maintenance
        J = Float64(0.), # maturity maintenance 
        C_W = (p.glb.C_W), # external stressor concentrations
        D_G = MVector{length(p.spc.k_D_G), Float64}(zeros(length(p.spc.k_D_G))), # scaled damage | growth efficiency
        D_M = MVector{length(p.spc.k_D_M), Float64}(zeros(length(p.spc.k_D_M))), # scaled damage | maintenance costs 
        D_A = MVector{length(p.spc.k_D_A), Float64}(zeros(length(p.spc.k_D_A))), # scaled damage | assimilation efficiency
        D_R = MVector{length(p.spc.k_D_R), Float64}(zeros(length(p.spc.k_D_R))), # scaled damage | reproduction efficiency
        D_h = MVector{length(p.spc.k_D_h), Float64}(zeros(length(p.spc.k_D_h))), # scaled damage | hazard rate
        y_G = Float64(1.), # relative response | growth efficiency
        y_M = Float64(1.), # relative response | maintenance costs 
        y_A = Float64(1.), # relative response | assimilation efficiency
        y_R = Float64(1.), # relative response | reproduction efficiency
        h_z = Float64(0.), # hazard rate | chemical stressors
        h_S = Float64(0.)  # hazard rate | starvation
    )
end



"""
    simulator(
        p::Ref{DEBParamCollection}; 
        alg = Tsit5(),
        saveat = 1, 
        reltol = 1e-6,
        kwargs...
        )::DataFrame

Run the DEBBase model from a `DEBParamCollection`  instance. <br>
These are the common parameters `p`. The agent-specific parameters `pagnt` are initialized by the simulator. <br>
Additional kwargs are passed on to `DifferentialEquations.solve()`.
"""
function simulator(
    p::DEBParamCollection; 
    system = DEB!,
    alg = Tsit5(),
    saveat = 1, 
    reltol = 1e-6,
    kwargs...
    )::DataFrame


    p.agn = AgentParams(p.spc) # initialize agent parameters incl. individual variability
    
    u = initialize_statevars(p)
    prob = ODEProblem(system, u, (0, p.glb.t_max), p) # define the problem
    sol = solve(prob, alg; saveat = saveat, reltol = reltol, kwargs...) # get solution to the IVP
    simout = sol_to_df(sol) # convert solution to dataframe
  
    return simout
end

"""
Run the DEBBase model from a reference to `DEBParamCollection`.
$(TYPEDSIGNATURES)
"""
function simulator(
    p::Ref{DEBParamCollection};
    saveat = 1,
    abstol = 1e-10, 
    reltol = 1e-10,
    kwargs...
    )

    return simulator(p.x; saveat = saveat, abstol = abstol, reltol = reltol, kwargs...) # run simulation
end


"""
    @replicates(simcall::Expr, nreps::Int64) 

Perform replicated runs of `simcall`, where `simcall` is a call to a simulator function. 

Example:

    spc = SpeciesParams(Z = Truncated(Normal(1, 0.1), 0, Inf)) # initialize default parameters with variable zoom factor
    yhat = @replicates DEBBase.simulator(DEBParamCollection(spc = spc))) 10 # execute replicated runs to simulator

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


"""
    sweep(simcall::Expr, component::AbstractParams, param::Symbol, range::Union{U,AbstractVector}) where {U <: UnitRange}

Perform a parameter sweep over a single parameter. 

- `simcall::Expr`: expression to evaluate a single iteration.
- `component::AbstractParams`: instance of a parameter struct which contains the parameter to be screened.
- `param::Symbol`: parameter to be screened.
- `range`: parameter range

----

Example:

    theta = DEBParamCollection() # use default parameters 
    theta.spc.Z = Truncated(Normal(1, 0.05), 0, Inf) # include agent variability

    @time yhat = sweep(
        :(@replicates simulator(theta) 10), # evaluate 10 replicates per iteration
        theta.spc, :Idot_max_rel, # screen parameter Idot_max_rel contained in theta.spc
        range(10, 25, length = 10) # evaluate 10 values between 10 and 25
        )
"""
function sweep(simcall::Expr, component::AbstractParams, param::Symbol, range::Union{U,AbstractVector}) where {U <: UnitRange}
    yhat = DataFrame()

    for val in range
        setproperty!(component, param, val)
        yhat_i = eval(simcall)
        yhat_i[!,param] .= val
        append!(yhat, yhat_i)
    end
    return yhat
end


"""
    @compose(derivs)
    
Compose a model system from a list of derivative functions. <br>
Each `deriv` is called as `deriv!(du, u, p..., t).`
"""
macro compose(derivs)
    quote
        function du!(du, u, p, t)
            for deriv! in $(esc(derivs))
                deriv!(du, u, p, t)
            end
        end
    end
end


