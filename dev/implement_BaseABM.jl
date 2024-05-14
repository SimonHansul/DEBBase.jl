@time "Loading packages" begin
    TAG = splitpath(@__FILE__)[end] |> x -> split(x, ".")[1] |> String    
    using Pkg; Pkg.activate("test")
    using Test
    using Plots, StatsPlots, Plots.Measures
    using SHUtils
    using DataFrames, DataFramesMeta
    using ProgressMeter
    using Distributions
    default(leg = false, lw = 1.5)

    using Revise
    using DEBBase
    using Parameters, DEBParamStructs
    using ComponentArrays
    include("../src/ModelFunctions.jl");
    include("../src/Simulators.jl");
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
## Model object
=#

"""
Definition of basic ABM object. <br>
Currently assumes that only a single species of type `AgentType` is simulated at a time.
"""
@with_kw mutable struct ABM <: AbstractABM
    p::AbstractParamCollection # parameter collection
    agents # agents
    t::Float64 # current simulation time
    dt::Float64 # timestep
    unique_id_count::Int64 # cumulative number of agents in the simulation
    u::ComponentVector # global state variables
    du::ComponentVector # global derivatives


    """
    Instantiate ABM from param collection `p`. 
    """
    function ABM(p::A; dt = 1/24) where A <: AbstractParamCollection
        abm = new() # initialize ABM object
        abm.p = p # store the parameters
        abm.t = 0. # initialize simulation time
        abm.dt = dt # timestep is given as keyword argument

        abm.u = init_substates_global(p) # initialize the global substates
        abm.du = similar(abm.u) # initialize the global derivatives
        init_agents!(abm) # initailize the population of agents

        return abm
    end
end


#=
## Agent Object
=#

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
    function BaseAgent(abm::AbstractABM)
        a = new()

        a.unique_id = abm.unique_id_count # identifier
        a.p = copy(p) # parameters
        a.p.agn = AgentParams(a.p.spc) # assign agent-level parameters (induces individual variability)
        initialize_statevars!(a, abm) # initialize agent-level state variables
        a.du = similar(a.u) # derivatives
        return a
    end
end


"""
    initialize_statevars!(agent::AbstractAgent)
Initialize agent-level state variables.
"""
function initialize_statevars!(agent::AbstractAgent, abm::AbstractABM)::Nothing

    agent.u = ComponentArray(
        glb = Ref(abm.u),
        agn = init_substates_agent(abm.p)
    )
    
    return nothing
end

"""
    initialize_statevars!(abm::AbstractABM)
Initialize ABM-level state variables.
"""
function initialize_statevars!(abm::AbstractABM)
    abm.u = ComponentArray(
        X_p = Float64(p.glb.Xdot_in), # initial resource abundance equal to influx rate
        C_W = (p.glb.C_W), # external stressor concentrations
    )
end

"""
    init_agents!(abm::AbstractABM)::nothing
Initialize the population of agents.
"""
function init_agents!(abm::AbstractABM)::nothing

    abm.agents = Vector{abm.p.glb.AgentType}(undef, p.glb.N0) # initialize a vector of agents with undefined values and defined length

    for i in 1:abm.p.glb.N0 # for the number of initial agents
        abm.agents[i] = BaseAgent(abm) # initialize an agent and add it to the vector of agents
    end

    return nothing
end

ABM(p)

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

    for (du_i, u_i, var) in zip(agent.du, agent.u, keys(agent.du))
        u_i = euler(du_i, u_i, abm.dt)
        setproperty!(agent, var, u_i)
    end
end

"""
    euler!(u::Float64, du::Float64, dt::Float64)
Apply Euler scheme.
"""
function euler(du::Float64, u::Float64,dt::Float64)
    return u + du * dt
end

"""
    euler!(u::Vector{Float64}, du::Vector{Float64}, dt::Float64)
Apply Euler scheme do a Vector of states and derivatives.
"""
function euler!(du::Vector{Float64}, u::Vector{Float64}, dt::Float64)
    return u .+ (du .* dt)
end


using ProfileView 

begin    
    p = DEBParamCollection()
    p.glb.N0 = 1_00
    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
    abm = ABM(p, dt = 1)

    aout = []

    @time while abm.t <= abm.p.glb.t_max
        for agent in abm.agents
            step!(agent, abm)

            #if (abm.t%1) == 0
            push!(aout, (t = abm.t, unique_id = agent.unique_id, u = agent.u))
            #end
        end

        abm.t += abm.dt # advance simulation time 
    end

    t = [a.t for a in aout]
    AgentID = [a.unique_id for a in aout]
    S = [a.u.S for a in aout]
    plot(t, S, group = AgentID)
    #plot(t, S, group = AgentID)
    #println(S)
end




#=
## Organization of state variables

within ABM constructor:

```
    abm.u_glb = ComponentVector(X_p = 0., C_W =  [0.])
    abm.du_glb = similar(u_glb)

    abm.u = ComponentVector(glb = Ref(abm.u_glb))
    abm.du = ComponentVector()
```


within agent constructor: 

```
    agent.u = ComponentVector(
        glb = Ref(abm.u_glb), 
        agn = ComponentVector(
            H = 0,
            ...
        )
        ...
    )

    agent.du = ComponentVector(
        glb = Ref(du_glb),
        agn = similar(agent.u.agn)

    )

```


=#