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
    using ComponentArrays
    using StaticArrays

    using Revise
    using DEBBase
    using Parameters, DEBParamStructs
    include("../src/ModelFunctions.jl");
    include("../src/Simulators.jl");
end

#=
## Implementation strategy 

- Ignore ODE simulation, everything is centered around ABM 
- Re-implement ODE later
    - Define ODESimulator, which initializes a single agent and uses an ODE solver
=#

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

    y!(du, u, p, t) # stressor responses

    Idot!(du, u, p, t) # ingestion flux
    Adot!(du, u, p, t) # assimilation flux
    Mdot!(du, u, p, t) # somatic maintenance flux
    Jdot!(du, u, p, t) # maturity maintenance flux

    Sdot!(du, u, p, t) # structure
    Hdot!(du, u, p, t) # maturity 
    H_bdot!(du, u, p, t) # estimate of maturity at birth
    Rdot!(du, u, p, t) # reproduction buffer
    
    X_pdot_out!(du, u, p, t) # resource abundance
    X_embdot!(du, u, p, t) # vitellus
    Ddot!(du, u, p, t) # damage
    C_Wdot!(du, u, p, t) # external stressor concentration

    for (du_i, u_i, var) in zip(agent.du.agn, agent.u.agn, keys(agent.du)) # for every state variable
        setproperty!(agent, var, euler(du_i, u_i, abm.dt)) # update agent substates
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
using DEBBase
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

using DEBBase

p = DEBParamCollection()
a = ABM(p)


u
a.du
a.agents[1].du.glb




