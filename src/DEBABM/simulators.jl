"""
    simulator(
        p::Union{NamedTuple,AbstractParamCollection}; 
        dt = 1/24, 
        saveat = 1,
        showinfo::Number = Inf
        )::DataFrame

Simulate the agent-based model from a parameter collection. 

```
p = ABMParamCollection()
sim = simulator(p)
```

args

- `p`: The parameter collection with defined global and species parameters.

kwargs

- `dt`: Length of a timestep in the model (unit according to chosen unit of rate parameters)
- `saveat`: Time interval at which to record output
- `showinfo`: Time interval at which to print an update. Nothing will be printed if `showinfo == Inf` (the default).

"""
function simulator(
    p::Union{NamedTuple,AbstractParamCollection}; 
    dt = 1/24, 
    saveat = 1,
    AgentParamType::DataType = AgentParams,
    showinfo::Number = Inf
    )::DataFrame

    showinfo < Inf ? @info("Running ABM simulation with t_max=$(p.glb.t_max)") : nothing
    
    m = ABM(p; dt = dt, saveat = saveat, AgentParamType = AgentParams)

    while !(m.t > m.p.glb.t_max)
        if showinfo < Inf && isapprox(m.t % showinfo, 0, atol = m.dt)
            @info("t=$(m.t)") : nothing
        end

        model_step!(m)
    end

    return agent_record_to_df(m)
end

"""
agent_record_to_df(
    m::AbstractDEBABM; 
    cols::Union{Symbol,Vector{Symbol}} = :all
    )::DataFrame

Convert agent record to Data Frame.

args 

- `m::AbstractDEBABM`: Model object containing an `agent_record` field as Vector of Component Vectors

kwargs

- `cols`: A Vector of Symbols indicating the 

"""
function agent_record_to_df(
    m::AbstractDEBABM; 
    )::DataFrame


    cols = vcat(
        [:t, :id],
        [keys(m.agent_record[1])...]
    )  |> unique

    
    hcat([map(x -> getproperty(x, y), m.agent_record) for y in cols]...) |> 
    x -> DataFrame(x, cols)
end
