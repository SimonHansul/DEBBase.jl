#=
# Base ABM implementation 

Implementation of the simplest possible generic version.<br>
Generic in terms of state variables and parameters, i.e. we don't need to make any chagnes to this code if we add or remove variables or parameters
=#

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
    using OrdinaryDiffEq
    using Agents
    using DEBBase
    using Parameters, DEBParamStructs
    using Revise
    include("../src/ModelFunctions.jl");
    include("../src/Simulators.jl");
    include("../src/IO.jl");
end


"""
    step!(agent::AbstractAgent, abm::DEBABM; odefuncs::Vector{Function}, rulefuncs::Vector{Function})

Definition of agent step. \n
This definition is generic, so that the function body does not have to be modified 
in order to simulate different models. \n 
An agent step is defined as the combination of an ODE system and rule-based functions. 
The definition of this system can be modified through the keyword argument `odefuncs` and `rulefuncs`.
`odefuncs` make up the part of the agent step which is represented as deterministic ordinary differential equations. 
`rulefuncs` make up the part of the agent step which is represented as rule-based functions. 
`odefuncs` and `rulefuncs` differ in their signature, but both are mutating:

- `odefunc(du, u, p, t)::Nothing`
- `rulefunc(agent, abm)::Nothing`


Args:

- `agent::AbstractAgent`: Any agent object
- `abm::DEBABM`: Any ABM object

Kwargs:

- `odefuncs::Vector{Function}`: Vector of ODE-based functions
- `rulefuncs::Vector{Function}`: Vector of rule-based functions
"""
function agent_step!(
    agent::DEBAgent, 
    abm::DEBABM, 
    odefuncs = [
        y!, 
        Idot!, 
        Adot!, 
        Mdot!, 
        Jdot!, 
        Sdot!, 
        Hdot!,
        H_bdot!,
        Rdot!,
        X_pdot_out!,
        X_embdot!,
        Ddot!
        ],
    rulefuncs = [
        reproduce!,
        die!
    ]
    )::Nothing

    du, u, p = agent.du, agent.u, agent.p # unpack agent substates, derivatives and parameters
    t = abm.t

    for func! in odefuncs # apply all ODE-based functions
        func!(du, u, p, t) 
    end

    for func! in rulefuncs # apply all rule-based functions
        func!(agent, abm)
    end

    map!(abm.euler, u.agn, du.agn, u.agn) # apply the euler scheme to agent substates

    return nothing
end


"""
    record!(agent::BaseAgent, abm::ABM)

Record agent data (only if `p.glb.recordagentvars == true`).
"""
function record!(agent::AbstractAgent, abm::DEBABM)::Nothing
    if abm.p.glb.recordagentvars && isapprox(abm.t%abm.saveat, 0, atol = abm.dt)
        push!(abm.aout, (t = abm.t, AgentID = agent.AgentID, u = copy(agent.u.agn)))
    end

    return nothing
end


"""
    step!(abm::ABM; odefuncs = [C_Wdot!, X_pdot_glb!], rulefuncs = [])::Nothing

Execute an ABM step. 

First, global functions are executed. 
As for the agent step, these are divided into ode-based functions and rule-based functions with the corresponding signatures: 

- `odefunc(du, u, p, t)::Nothing`
- `rulefunc(abm)::Nothing`

Secondly, the agent steps are executed.

"""
function model_step!(
    abm::DEBABM; 
    odefuncs = [C_Wdot_const!, X_pdot_chemstat!], 
    rulefuncs = [])::Nothing

    du, u, p = abm.du, abm.u, abm.p
    t = abm.t

    for func! in odefuncs # execute global ode-based functions
        func!(du, u, p, t)
    end

    for func! in rulefuncs # execute global rule-based functions
        func!(abm)
    end

    for a in filter(a -> !a.dead, abm.agents) # for every living agent
        agent_step!(a, abm) # execute agent step
        record!(a, abm) # record agent data
    end

    abm.t += abm.dt

    return nothing
end


p = DEBParamCollection()



    StandardABM(
    DEBAgent; 
    agent_step! = agent_step!, 
    model_step! = model_step!,
    properties = DEBABM
    )




