function init_substates_agent()
    return ComponentArray( # initial states
        X_emb = Float64(p.agn.X_emb_int), # initial mass of vitellus
        S = Float64(p.agn.X_emb_int * 0.001), # initial structure is a small fraction of initial reserve // mass of vitellus
        H = Float64(0), # maturity
        H_b = Float64(0), # maturity at birth (will be derived from model output)
        R = Float64(0), # reproduction buffer
        I_emb = Float64(0), # ingestion from vitellus
        I_p = Float64(0), # ingestion from external food resource
        I = Float64(0), # total ingestion
        A = Float64(0), # assimilation
        M = Float64(0), # somatic maintenance
        J = Float64(0), # maturity maintenance 

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
    initialize_statevars(p::AbstractParamCollection, pagnt::ComponentVector{Float64})::ComponentArray 
For initialization of ODE simulator, initialize the component vector of state variables, `u`, based on common parameters `p`.
"""
function initialize_statevars(p::AbstractParamCollection)::ComponentArray 
    glb = ComponentArray(
        x = ComponentArray( # we add the 'x' to mimic the ABM syntax where agn.u.glb is a reference to the actual glbinstance
            X_p = Float64(p.glb.Xdot_in), # initial resource abundance equal to influx rate
            C_W = (p.glb.C_W) # external stressor concentrations
        )
    )

    agn = init_substates_agent()

    return ComponentArray(
        glb = glb,
        agn = agn
    )
end

@enum ReturnType dataframe odesolution matrix # possible return types

function abstractsimulator(
    p::AbstractParamCollection,
    model, 
    AgentParamType::DataType;
    alg = Tsit5(),
    returntype::ReturnType = dataframe,
    kwargs...
    )::Union{DataFrame,ODESolution}

    p.agn = AgentParamType(p.spc) # initialize agent parameters incl. individual variability
    
    u = initialize_statevars(p)
    prob = ODEProblem(model, u, (0, p.glb.t_max), p) # define the problem
    sol = solve(prob, alg; kwargs...) # get solution to the IVP

    if returntype == dataframe
        simout = sol_to_df(sol) # convert solution to dataframe
        return simout
    end

    if returntype == odesolution
        return sol
    end

    if returntype == matrix
        simout = sol_to_mat(sol)
        return simout
    end
end

"""
    simulator(
        p::Ref{DEBParamCollection}; 
        alg = Tsit5(),
        saveat = 1, 
        reltol = 1e-6,
        kwargs...
        )::DataFrame

Run the DEBBase model from a `DEBParamCollection` instance. <br>
"""
function simulator(
    p::DEBParamCollection; 
    model = DEB!,
    AgentParamType::DataType = AgentParams,
    kwargs...
    )

    abstractsimulator(
        p, 
        model,
        AgentParamType;     
        alg = Tsit5(),
        saveat = 1, 
        reltol = 1e-6,
        kwargs...
        )
end

"""
    simulator(
        p::Ref{DEBParamCollection};
        saveat = 1,
        kwargs...
        )
Run the DEBBase model from a reference to `DEBParamCollection`.
"""
function simulator(
    p::Ref{DEBParamCollection};
    saveat = 1,
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

    yhat = sweep(
        :(@replicates simulator(theta) 10), # evaluate 10 replicates per iteration
        theta.spc, :Idot_max_rel, # screen parameter Idot_max_rel contained in theta.spc
        range(10, 25, length = 10) # evaluate 10 values between 10 and 25
        )
"""
#macro sweep(
#    simcall::Expr, 
#    component::A,
#    param::Symbol, 
#    range::Union{U,AbstractVector}) where {U <: UnitRange, A <: AbstractParams}
#    
#    quote
#        yhat = DataFrame()
#
#        for val in $range
#            setproperty!($component, $param, $val)
#            yhat_i = $(esc(simcall))
#            yhat_i[!,$param] .= val
#            append!(yhat, yhat_i)
#        end
#        yhat
#    end
#end


#"""
# NOTE : This macro is currently outcommented because it turned out that models defined with @compose execute very slowly.
# Also, the usefulness of @compose has to be critically evaluated.
#
#    @compose(derivs)
#    
#Compose a model system from a list of derivative functions. <br>
#Each `deriv` is called as `deriv!(du, u, p..., t).`
#"""
#macro compose(derivs)
#    quote
#        function du!(du, u, p, t)
#            for deriv! in $(esc(derivs))
#                deriv!(du, u, p, t)
#            end
#        end
#    end
#end


