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
    using ComponentArrays
    include("../src/ModelFunctions.jl")
end

#=
We first (again) test the @replicates macro, 
which also tests the initalization of agent params and induction of individual variability.
=#

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

using DEBBase

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
        a.p = copy(p) # parameters
        a.p.agn = AgentParams(a.p.spc) # assign agent-level parameters (induces individual variability)
        a.u = DEBBase.initialize_statevars(a.p) # state variables
        a.du = similar(a.u) # derivatives
        return a
    end
end

"""
    init_pop(p::DEBParamCollection, N0::Int64; AgentType = BaseAgent)
Initialize a 
"""
function init_pop(p::DEBParamCollection; AgentType = BaseAgent)
    agents = Vector{AgentType}(undef, p.glb.N0)

    for i in 1:p.glb.N0
        a = BaseAgent(p, i)
        agents[i] = a
    end

    return agents
end

abstract type AbstractABM end

"""
Definition of basic ABM object. <br>
Currently assumes that only a single species of type `AgentType` is simulated at a time.
"""
@with_kw mutable struct ABM <: AbstractABM
    p::AbstractParamCollection # Parameter collection
    agents # Agents
    t::Float64 # current simulation time
    dt::Float64 # timestep

    """
    Instantiate ABM from param collection `p`. 
    """
    function ABM(theta::A; AgentType = BaseAgent, dt = 1/24) where A <: AbstractParamCollection
        abm = new() # initialize ABM object
        abm.p= p # 
        abm.t = 0.
        abm.dt = dt

        abm.agents = init_pop(p)

        return abm
    end
end

p = DEBParamCollection()
p.glb.N0 = 10
p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
@time abm = ABM(p)

#=
Agent step
=#

"""
    step!(agent::AbstractAgent, abm::AbstractABM)
Definition of agent step. 
"""
function step!(agent::AbstractAgent, abm::AbstractABM)
    du, u, p = agent.du, agent.u, agent.p 
    t = abm.t

    #### stressor responses
    y!(du, u, p, t)

    #### auxiliary state variables (record cumulative values)
    Idot!(du, u, p, t)
    Adot!(du, u, p, t) 
    Mdot!(du, u, p, t) 
    Jdot!(du, u, p, t)

    #### major state variables
    Sdot!(du, u, p, t) # structure
    Hdot!(du, u, p, t) # maturity 
    H_bdot!(du, u, p, t) # estimate of maturity at birth
    Rdot!(du, u, p, t) # reproduction buffer
    
    X_pdot_out!(du, u, p, t) # resource abundance
    X_embdot!(du, u, p, t) # vitellus
    Ddot!(du, u, p, t) # damage
    C_Wdot!(du, u, p, t) # external stressor concentration

    update!(agent, abm)
end

function euler!(u::Float64, du::Float64, dt::Float64)
    return u + du * dt
end

function euler!(u::Vector{Float64}, du::Vector{Float64}, dt::Float64)
    return euler!.(u, du, dt)
end

begin
    function update!(agent::AbstractAgent, abm::AbstractABM)
        for (du_i, u_i) in zip(agent.du, agent.u)
            u_i = euler!(u_i, du_i, abm.dt)
        end
    end
    
    abm = ABM(p)

    aout = []

    while abm.t <= abm.p.glb.t_max
        for agent in abm.agents
            step!(agent, abm)
            push!(aout, (t = abm.t, unique_id = agent.unique_id, u = agent.u))
        end

        abm.t += abm.dt
    end

    # FIXME: u is not updated...

    t = [x.t for x in aout]
    unique_id = [x.unique_id for x in aout]
    S = [x.u.S for x in aout]
    
    plot(t, S, group = unique_id)
end




p.glb