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
    using Parameters, DEBParamStructs
    using Random
    using Agents

    using Revise
    using DEBBase
    
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
    step!(agent::Agents.AbstractAgent, abm::AbstractABM; odefuncs::Vector{Function}, rulefuncs::Vector{Function})

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

- `agent::Agents.AbstractAgent`: Any agent object
- `abm::AbstractABM`: Any ABM object

Kwargs:

- `odefuncs::Vector{Function}`: Vector of ODE-based functions
- `rulefuncs::Vector{Function}`: Vector of rule-based functions
"""
function step!(
    agent::BaseAgent, abm::ABM, 
    odefuncs = [
        y_z!, 
        h_S!,
        Idot!, 
        Adot!, 
        Mdot!, 
        Jdot!, 
        Sdot!, 
        Hdot!,
        H_bdot!,
        Rdot!,
        Ddot!,
        age!
        ],
    rulefuncs = [
        reproduce_opportunistic!,
        die!
    ])

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

Record agent output (only if `p.glb.recordagentvars == true`).
"""
function record!(agent::BaseAgent, abm::ABM)
    if abm.p.glb.recordagentvars && isapprox(abm.t%abm.saveat, 0, atol = abm.dt)
        push!(abm.aout, (t = abm.t, AgentID = agent.AgentID, u = copy(agent.u.agn)))
    end
end

"""
    record!(abm::ABM)

Record model-level output.
"""
function record!(abm::ABM)
    if isapprox(abm.t%abm.saveat, 0, atol = abm.dt)
        push!(abm.mout, (t = abm.t, u = copy(abm.u)))
    end
end

function N_tot!(abm::AbstractABM)::Nothing

    abm.u.N_tot = length(abm.agents)

    return nothing
end

"""
    step!(
        abm::ABM; 
        sysfuncs = [C_Wdot_const!, X_pdot_chemstat!],
        rulefuncs = [N_tot!]
        )

Execution of a generic ABM step, following the schedule: 

1. Shuffle agents
2. Calculate global derivatives
3. Calculate agent derivatives
4. Update agent states
5. Record agent states
6. Update global states
7. Record global states

It is important that the agent steps occur between calculation of the global derivatives and 
updating the global states, because global derivatives which may be influenced by agent derivatives 
are initialized during calculation of the global states and mutated by the agents. 
Changing this order will lead to erroneous computation of the interaction between agents and the environment.
"""
function step!(
    abm::ABM; 
    sysfuncs = [C_Wdot_const!, X_pdot_chemstat!],
    rulefuncs = [N_tot!]
    )
    du, u, p = abm.du, abm.u, abm.p
    t = abm.t

    shuffle!(abm.agents)

    for func! in sysfuncs # execute global functions
        func!(du, u, p, t)
    end
    
    for a in abm.agents # for every agent
        step!(a, abm) # execute agent steps
        record!(a, abm) # record agent data
    end
    
    map!(abm.euler, u, du, u) # apply the euler scheme to global states
    record!(abm)
    filter!(a -> a.u.agn.dead == false, abm.agents) # remove dead agents

    abm.t += abm.dt
end

function run!(abm::AbstractABM)
    while abm.t <= abm.p.glb.t_max
        step!(abm)
    end
end

#=
##FIXME: food density is not updated correctly..

We shouldn't need reference to u.glb / du.glb as long as nothing is copied explicitly...

- Replacing Ref(glb) with glb 
- Replacing glb.x with glb

=#

begin    
    using DEBBase, DoseResponse
    using DEBFigures
    include("../src/ModelFunctions.jl")
    
    p = DEBParamCollection()
    p.glb.N0 = 1

    p.glb.k_V = 0.5
    p.glb.Xdot_in = 10.
    p.glb.t_max = 20.
    p.spc.e_S = 0.5
    p.spc.b_S = 100.
    p.spc.eta_AR = 0.

    function simulator(p::AmphiDEBParamCollection)

    abm = ABM(p)
    @time run!(abm)
    
    aout_df = DataFrame(hcat([[x.t, x.AgentID] for x in abm.aout]...)', [:t, :AgentID]) |> 
    x-> hcat(x, DataFrame(hcat([a.u for a in abm.aout]...)', extract_colnames(abm.aout[1].u)))
    
    mout_df = DataFrame(hcat([[x.t] for x in abm.mout]...)', [:t]) |> 
    x-> hcat(x, DataFrame(hcat([m.u for m in abm.mout]...)', extract_colnames(abm.mout[1].u)))

    pa = @df aout_df plot(
        groupedlineplot(:t, :S, :cohort)
    )
    
    pN = combine(groupby(aout_df, :t), :AgentID => x -> length(unique(x))) |> 
    x -> @df x plot(:t, :AgentID_function)
    
    pm = @df mout_df plot(:t, :X_p)

    plot(pa, pm, pN)
end   

#=
 
=#


@df aout_df plot(:t, :S ./ :S_0, group = :AgentID)


@agent struct DEBAgent(NoSpaceAgent)
    p::AbstractParamCollection
    u::ComponentVector
    du::ComponentVector
    AgentID::Int64

    """
        BaseAgent(p::AbstractParamCollection, AgentID::Int64)
    
    Initialize a base agent from parameters. 
    """
    function BaseAgent(abm::ABM)
        a = new()

        a.p = copy(abm.p) # parameters; TODO: #29 avoid copying everything
        a.p.agn = AgentParams(a.p.spc) # assign agent-level parameters (induces individual variability)
        initialize_statevars!(a, abm) # initialize agent-level state variables
        a.du = similar(a.u) # initialize agent-level derivatives
        a.du.glb = abm.du # reference to global derivatives
        a.du.agn = copy(a.u.agn) # derivatives of agent substates have same shape as states
        a.AgentID = abm.AgentID_count # assign AgentID
        abm.AgentID_count += 1 # increment AgentID counter
        
        return a
    end
end