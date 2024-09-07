#"""
#Composes an ODE system from a vector of global and species-specific derivatives, respectively. 
#This brings two advantages: 
#
#- Function definitions can be recycled in different simulation contexts (e.g. ODE-based and IBM)
#- Function vectors can be combined to facilitate modular modelling: Coupling two modules is equivalent to concatenating the function vectors and adding an adapter
#
#A potential pitfall of this approach is to split up the model into too many small functions, making the model less readable. 
#
#Users generally don't need to call with `compositemodel`, this function is used internally in `simulator`.
#"""
#function compositemodel!(du, u, p, t)::Nothing
#
#    for func! in p.glb.odefuncs # calculate the global derivatives
#        func!(du, u, p, t) 
#    end
#
#    for func! in p.spc.odefuncs # calculate the species-specific derivatives
#        func!(du, u, p, t)
#    end
#
#    return nothing
#end


"""
We use callbacks to keep track of the current life stage, avoiding if/else statements in 
the derivatives. 
Each life stage is associated with an identically named auxiliary variable, 
taking logical values to indicate whether this is the current life stage.
"""
function lifestage_callbacks()

    function condition_embryo(u, t, integrator) 
        u.X_emb < 0
    end

    function effect_embryo!(integrator) 
        integrator.u.embryo = 1.
        integrator.u.juvenile = 0.
        integrator.u.adult = 0.
    end
    cb_embryo = DiscreteCallback(condition_embryo, effect_embryo!)

    condition_juvenile(u, t, integrator) = (u.X_emb <= 0) && (u.H < u.H_p)
    function effect_juvenile!(integrator) 
        integrator.u.embryo = 0.
        integrator.u.juvenile = 1.
        integrator.u.adult = 0.
    end
    cb_juvenile = DiscreteCallback(condition_juvenile, effect_juvenile!)

    condition_adult(u, t, integrator) = u.H >= u.H_p
    function effect_adult!(integrator) 
        integrator.u.embryo = 0.
        integrator.u.juvenile = 0.
        integrator.u.adult = 1.
    end
    cb_adult = DiscreteCallback(condition_adult, effect_adult!)

    return CallbackSet(cb_embryo, cb_juvenile, cb_adult)
end

"""
simulator(
    p::Union{AbstractParamCollection,NamedTuple}; 
    model = DEBODE!,
    alg = Tsit5(),
    saveat = 1,
    reltol = 1e-6,
    AgentParamType::DataType = ODEAgentParams,
    kwargs...)::DataFrame

Run an ODE-based model. 

**args**:

- `theta::Union{AbstractParamCollection,NamedTuple}`: A parameter collection with defined global parameters (`<: AbstractGlobalParams`) and species parameters (`<: AbstractSpeciesParams`).

**kwargs**:

- `model = DEBODE!`: Definition of the derivatives. A function form `du!(du, u, p, t)::Nothing`. See definition of `DEBODE!` in `derivatives.jl` for an example, or see docs and examples of [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) for more details.
- `alg = Tsit5()`: Algorithm to be used by `solve` function 
- `saveat = 1`: When or how often to save ODE solutions to output
- `reltol = 1e-6`: Relative tolerance of ODE solutions
- `AgentParamType::DataType = ODEAgentParams`: The data type that stores those parameters that are affected by individual variability. There has to be a corresponding constructor so that `theta.agn = AgentParamType(theta.spc)` works. 
- `kwargs...`: Additional keyword argument are passed on to `OrdinaryDiffEq.solve`

**Example**: 

```Julia
yhat = simulator(DEBParamCollection())
```

"""
function simulator(
    theta::Union{AbstractParamCollection,NamedTuple}; 
    alg = Tsit5(),
    saveat = 1,
    reltol = 1e-6,
    model = DEBBase!,
    AgentParamType::DataType = ODEAgentParams,
    kwargs...
    )::DataFrame

    theta.agn = AgentParamType(theta.spc) # initialize agent parameters incl. individual variability
    callbacks = lifestage_callbacks()

    u = initialize_statevars(theta)
    prob = ODEProblem(model, u, (0, theta.glb.t_max), theta) # define the problem
    sol = solve(prob, alg; callback = callbacks, saveat = saveat, reltol = reltol, kwargs...) # get solution to the IVP
    simout = sol_to_df(sol) # convert solution to dataframe

    return simout
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