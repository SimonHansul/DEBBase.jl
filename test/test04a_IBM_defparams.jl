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
    function BaseAgent(abm::AbstractABM)
        a = new()

        a.unique_id = abm.unique_id_count # identifier
        a.p = copy(p) # parameters
        a.p.agn = AgentParams(a.p.spc) # assign agent-level parameters (induces individual variability)
        initialize_statevars!(a) # initialize agent-level state variables
        a.du = similar(a.u) # derivatives
        return a
    end
end


"""
    initialize_statevars!(agent::AbstractAgent)
Initialize agent-level state variables.
"""
function initialize_statevars!(agent::AbstractAgent)
    agent.u = ComponentArray( # initial states
        #X_p = Float64(p.glb.Xdot_in), # initial resource abundance equal to influx rate
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
        #C_W = (p.glb.C_W), # external stressor concentrations
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
    p::AbstractParamCollection # parameter collection
    agents # agents
    t::Float64 # current simulation time
    dt::Float64 # timestep
    unique_id_count::Int64 # cumulative number of agents in the simulation

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