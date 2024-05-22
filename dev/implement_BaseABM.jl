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

    using DEBBase
    using Parameters, DEBParamStructs
    using Revise
    include("../src/ModelFunctions.jl");
    include("../src/Simulators.jl");
    include("../src/IO.jl");
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
    step!(agent::AbstractAgent, abm::AbstractABM; odefuncs::Vector{Function}, rulefuncs::Vector{Function})

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
- `abm::AbstractABM`: Any ABM object

Kwargs:

- `odefuncs::Vector{Function}`: Vector of ODE-based functions
- `rulefuncs::Vector{Function}`: Vector of rule-based functions
"""
function step!(
    agent::BaseAgent, abm::ABM, 
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
    )
    du, u, p = agent.du, agent.u, agent.p # unpack agent substates, derivatives and parameters
    t = abm.t

    for func! in odefuncs # apply all ODE-based functions
        func!(du, u, p, t) 
    end

    for func! in rulefuncs # apply all rule-based functions
        func!(agent, abm)
    end

    map!(abm.euler, u.agn, du.agn, u.agn) # apply the euler scheme to agent substates
end

"""
    record!(agent::BaseAgent, abm::ABM)

Record agent data (only if `p.glb.recordagentvars == true`).
"""
function record!(agent::BaseAgent, abm::ABM)
    if abm.p.glb.recordagentvars && isapprox(abm.t%abm.saveat, 0, atol = abm.dt)
        push!(abm.aout, (t = abm.t, AgentID = agent.AgentID, u = copy(agent.u.agn)))
    end
end

function step!(abm::ABM; sysfuncs = [C_Wdot!, X_pdot_glb!])
    du, u, p = abm.du, abm.u, abm.p
    t = abm.t

    for func! in sysfuncs # execute global functions
        func!(du, u, p, t)
    end

    for a in abm.agents 
        step!(a, abm) # execute agent steps
        record!(a, abm) # record agent data
    end
end


using DEBBase
using ProfileView
using BenchmarkTools

begin

    #=
    ## Test: Simulate agents independently
    =#

    using ProfileView
    using DEBBase

    p = DEBParamCollection()
    abm = ABM(p)

    @time begin    
        p = DEBParamCollection()
        p.glb.N0 = 100
        p.glb.t_max = 56.
        p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
        abm = ABM(p, dt = 1/24)


        aout = []

        while abm.t <= abm.p.glb.t_max
            for a in abm.agents
                step!(a, abm)

                if abm.p.glb.recordagentvars && isapprox(abm.t%abm.saveat, 0, atol = abm.dt)
                    push!(aout, (t = abm.t, AgentID = a.AgentID, u = copy(a.u.agn)))
                end
            end
            
            abm.t += abm.dt # advance simulation time 
        end
    end
end

aout_df = DataFrame(hcat([[x.t, x.AgentID] for x in aout]...)', [:t, :AgentID]) |> 
x-> hcat(x, DataFrame(hcat([a.u for a in aout]...)', extract_colnames(aout[1].u)))


@df aout_df plot(
