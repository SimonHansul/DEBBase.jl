@time "Loading packages" begin
    using Pkg; Pkg.activate("test")
    using Test
    using Plots, StatsPlots, Plots.Measures
    using SHUtils
    using DataFrames, DataFramesMeta
    using ProgressMeter
    using Distributions
    default(leg = false, lw = 1.5)
    TAG = splitpath(@__FILE__)[end] |> x -> split(x, ".")[1] |> String    

    using Revise
    using DEBBase
    using Parameters, DEBParamStructs
end

@test begin
    plt = plot(
        layout = (1,2), 
        leftmargin = 2.5mm, bottommargin = 2.5mm, 
        xlabelfontsize = 10, 
        size = (600,350),
        xlabel = "Time since fertilization (d)",
        ylabel = ["Structure" "Reproduction buffer"]
        )

    yhat_fix = DEBBase.simulator(DEBParamCollection())
    @df yhat_fix plot!(plt, :t, :S, subplot = 1, color = :black, lw = 2)
    @df yhat_fix plot!(plt, :t, :R, subplot = 2, color = :black, lw = 2)
    
    hyperZ = Truncated(Normal(1., 0.1), 0, Inf)
    yhat_var = @replicates DEBBase.simulator(DEBParamCollection(spc = SpeciesParams(Z = hyperZ))) 10
    @df yhat_var plot!(plt, :t, :S, group = :replicate, alpha = .25, subplot = 1, c = :viridis)
    @df yhat_var plot!(plt, :t, :R, group = :replicate, alpha = .25, subplot = 2, c = :viridis) 
    display(plt)
    true
end

#=
## Implementing an Agent object. 

What does an agent need?

- A reference to common paramters `pcmn`
- Agen-specific parameters `pown`
- State variables `u`
- Derivatives `du`

- A function to initialize the agent
=#

using ComponentArrays
abstract type AbstractAgent end

"""
DEBBase Agent. <br>
Each agent owns a reference to its associated parameter collection.
"""
mutable struct BaseAgent <: AbstractAgent
    p::AbstractParamCollection
    u::ComponentVector
    du::ComponentVector
    unique_id::Int64

    """
        BaseAgent(p::AbstractParamCollection, unique_id::Int64)
    
    Initialize a base agent from parameters. 
    """
    function BaseAgent(p::AbstractParamCollection, unique_id::Int64)
        a = new()

        a.unique_id = unique_id # identifier
        a.p = p # parameters
        a.p.agn = AgentParams(a.p.spc) # agent-level parameters (induces individual variability)
        a.u = DEBBase.initialize_statevars(a.p) # state variables
        a.du = similar(a.u) # derivatives
        return a
    end
end

p = DEBParamCollection()
p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
a = BaseAgent(p, 1)    

"""
    init_pop(p::DEBParamCollection, N0::Int64; AgentType = BaseAgent)
Initialize a 
"""
function init_pop(p::DEBParamCollection, N0::Int64; AgentType = BaseAgent)
    agents = Vector{AgentType}(undef, N0)

    for i in 1:N0
        agents[i] = BaseAgent(p, i)
    end

    return agents
end

agents = init_pop(p, 10)


abstract type AbstractABM end
mutable struct BaseABM <: AbstractABM
    X_p::Float64 # food concentration
    unique_id_count::Int64 # cumulative unique ids
    agents::Vector{AbstractAgent}

    function BaseABM(p::AbstractParamCollection)
    end
end







#=
## Implementing a model object 


=#



θ = DEBParamCollection()
θ_ref = Ref(θ)

agents = Vector{AbstractAgent}(undef, 10)
agents[1] = BaseAgent(θ_ref)

"""
Definition of basic ABM object. <br>
Currently assumes that only a single species of type `AgentType` is simulated at a time.
"""
@with_kw mutable struct ABM <: AbstractABM
    theta::AbstractParamCollection # Parameter collection
    agents # Agents
    t::Float64

    """
    Instantiate ABM from param collection `theta`. 
    """
    function ABM(theta::A; AgentType = BaseAgent, N_max = Int(1e4)) where A <: AbstractParamCollection
        abm = new() # initialize ABM object
        abm.theta = theta # 
        abm.t = 0.

        abm.agents = Vector{AbstractAgent}(undef, N_max)
        
        for i in eachindex(abm.agents)
            abm.agents[i] = AgentType(Ref(theta))
            if i <= theta.glb.N0
                activate!(abm.agents[i])
            end
        end

        return abm
    end
end


abm = ABM(θ)
println(abm)