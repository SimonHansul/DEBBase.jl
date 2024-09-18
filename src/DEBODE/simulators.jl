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

    condition_juvenile(u, t, integrator) = u.X_emb # transition to juvenile when X_emb hits 0
    function effect_juvenile!(integrator) 
        integrator.u.embryo = 0.
        integrator.u.juvenile = 1.
        integrator.u.adult = 0.
    end
    cb_juvenile = ContinuousCallback(condition_juvenile, effect_juvenile!)

    condition_adult(u, t, integrator) = u.H_p - u.H # condition to adult when H reaches H_p
    function effect_adult!(integrator) 
        integrator.u.embryo = 0.
        integrator.u.juvenile = 0.
        integrator.u.adult = 1.
    end
    cb_adult = ContinuousCallback(condition_adult, effect_adult!)

    return CallbackSet(cb_juvenile, cb_adult)
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
sim = simulator(DEBParamCollection())
```

"""
function simulator(
    params::Union{AbstractParamCollection,NamedTuple}; 
    alg = Tsit5(),
    saveat = 1,
    reltol = 1e-6,
    model = DEBODE_IA!,
    AgentParamType::DataType = ODEAgentParams,
    kwargs...
    )::DataFrame

    params.agn = AgentParamType(params.spc) # initialize agent parameters incl. individual variability
    callbacks = lifestage_callbacks()

    u = initialize_statevars(params)
    prob = ODEProblem(model, u, (0, params.glb.t_max), params) # define the problem
    sol = solve(prob, alg; callback = callbacks, saveat = saveat, reltol = reltol, kwargs...) # get solution to the IVP
    simout = sol_to_df(sol) # convert solution to dataframe

    return simout
end


"""
    @replicates(simcall::Expr, nreps::Int64) 

Perform replicated runs of `simcall`, where `simcall` is a call to a simulator function. 

Example:

```Julia
    spc = SpeciesParams(Z = Truncated(Normal(1, 0.1), 0, Inf)) # initialize default parameters with variable zoom factor
    sim = @replicates DEBBase.simulator(DEBParamCollection(spc = spc))) 10 # execute replicated runs to simulator
```

In this case, `sim` will contain the output of 10 replicated simulations. For each replicate, the zoom factor is sampled from a truncated Normal distribution. 
`sim` contains an additional column `replicate`.
"""
macro replicates(simcall::Expr, nreps::Int64)
    quote
        sim = DataFrame()

        for replicate in 1:$nreps
            sim_i = $(esc(simcall))
            sim_i[!,:replicate] .= replicate
            append!(sim, sim_i)
        end
        sim
    end
end


"""
    replicates(simulator::Function, params::Union{NamedTuple,AbstractParamCollection}, nreps::Int64; kwargs...)

Perform replicated runs of `simulator` with parameters `params` (`simulator(params)` has to be a valid function call). 
Analogous to `@replicates`, but a bit more flexible.
"""
function replicates(simulator::Function, params::Union{NamedTuple,AbstractParamCollection}, nreps::Int64; kwargs...)
    sim = DataFrame()

    for replicate in 1:nreps
        sim_i = simulator(params; kwargs...)
        sim_i[!,:replicate] .= replicate
        append!(sim, sim_i)
    end
    
    sim
end


"""
    exposure(simcall::Expr, C_Wvec::Vector{Float64}; kwargs...)

Simulate exposure to a single stressor over a Vector of constant exposure concentrations `C_Wvec`. 

"""
function exposure(simulator::Function, params::Union{AbstractParamCollection,NamedTuple}, C_Wvec::Vector{Float64})
    
    let C_W_int = params.glb.C_W # we will modify this value and then reset to the initial value
        sim = DataFrame()

        for C_W in C_Wvec
            params.glb.C_W[1] = C_W
            sim_i = simulator(params)
            append!(sim, sim_i)
        end

        params.glb.C_W = C_W_int 

        return sim
    end
end
