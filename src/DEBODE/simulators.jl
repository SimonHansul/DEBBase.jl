"""
Compose an ODE system from a vector of global and species-specific derivatives, respectively.
"""
function compositemodel!(du, u, p, t)::Nothing

    for func! in p.glb.odefuncs # calculate the global derivatives
        func!(du, u, p, t) 
    end

    for func! in p.spc.odefuncs # calculate the species-specific derivatives
        func!(du, u, p, t)
    end

    return nothing
end


"""
simulator(
    p::AbstractParamCollection; 
    model = DEBODE!,
    alg = Tsit5(),
    saveat = 1,
    reltol = 1e-6,
    AgentParamType::DataType = ODEAgentParams,
    kwargs...)::DataFrame

Run an ODE-based model. 

**args**:

- `theta::AbstractParamCollection`: A parameter collection with defined global parameters (`<: AbstractGlobalParams`) and species parameters (`<: AbstractSpeciesParams`).

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
    theta::AbstractParamCollection; 
    model = DEBODE!,
    alg = Tsit5(),
    saveat = 1,
    reltol = 1e-6,
    AgentParamType::DataType = ODEAgentParams,
    kwargs...
    )::DataFrame

    theta.agn = AgentParamType(theta.spc) # initialize agent parameters incl. individual variability

    u = initialize_statevars(theta)
    prob = ODEProblem(model, u, (0, theta.glb.t_max), theta) # define the problem
    sol = solve(prob, alg; saveat = saveat, reltol = reltol, kwargs...) # get solution to the IVP
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