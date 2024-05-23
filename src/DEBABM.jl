

"""
step!(agent::AbstractAgent, abm::AbstractABM; odefuncs::Vector{Function}, rulefuncs::Vector{Function})

Definition of agent step. \n
This definition is generic, so that the function body does not have to be modified 
in order to simulate different models.
"""
function step!(agent::BaseAgent, abm::ABM)
    du, u, p = agent.du, agent.u, agent.p # unpack agent substates, derivatives and parameters
    t = abm.t

    for func! in agent.odefuncs # apply all ODE-based functions
        func!(du, u, p, t) 
    end

    for func! in agent.rulefuncs # apply all rule-based functions
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
are initialized during calculation of the global states and mutated by the  
Changing this order will lead to erroneous computation of the interaction between agents and the environment.
"""
function step!(abm::ABM)

    du, u, p = abm.du, abm.u, abm.p
    t = abm.t

    shuffle!(abm.agents)

    for func! in abm.odefuncs # execute global ODE-based functions
        func!(du, u, p, t)
    end

    for func! in abm.rulefuncs # execute global rule-based functions
        func!(abm)
    end

    for a in abm.agents # for every agent
        step!(a, abm) # execute agent steps
        record!(a, abm) # record agent data
    end

    map!(abm.euler, u, du, u) # apply the euler scheme to global states
    record!(abm) # record global states
    filter!(a -> a.u.agn.dead == false, abm.agents) # remove dead agents

    abm.t += abm.dt
end

function run!(abm::AbstractABM)
    while abm.t <= abm.p.glb.t_max
        step!(abm)
    end
end

"""
Definition of basic ABM object. <br>
Currently assumes that only a single species of type `AgentType` is simulated at a time.
"""
mutable struct ABM <: AbstractABM

    odefuncs::Vector{Function} # ODE-based model step functions
    rulefuncs::Vector{Function} # rule-based model step functions

    p::AbstractParamCollection # parameter collection
    t::Float64 # current simulation time
    dt::Float64 # timestep
    euler::Function # definition of the euler function for the given dt
    saveat::Float64
    
    u::ComponentVector # global state variables
    du::ComponentVector # global derivatives

    AgentID_count::Int64 # cumulative number of agents in the simulation
    agents::AbstractVector # population of agents
    aout::AbstractVector # agent output
    mout::AbstractVector # model output

    """
    Instantiate ABM from param collection `p`. 
    """
    function ABM(p::A; dt = 1/24, saveat = 1) where A <: AbstractParamCollection
        abm = new() # initialize ABM object

        abm.odefuncs = Function[C_Wdot_const!, X_pdot_chemstat!] 
        abm.rulefuncs = Function[N_tot!]

        abm.p = p # store the parameters
        abm.t = 0. # initialize simulation time
        abm.dt = dt # timestep is given as keyword argument
        abm.euler = defeuler(dt) # Euler function for the given dt
        abm.saveat = saveat
        
        abm.u = init_substates_global(p) # initialize the global substates
        abm.du = similar(abm.u) # initialize the global derivatives
        abm.du.X_p = 0.
        abm.du.C_W .= 0.
        
        abm.AgentID_count = 0 # set agent count
        initialize_agents!(abm) # initailize the population of agents
        abm.aout = [] # initialize agent output
        abm.mout = [] # initialize model output

        run!(abm)

        aout = DataFrame(hcat([[x.t, x.AgentID] for x in abm.aout]...)', [:t, :AgentID]) |> 
        x-> hcat(x, DataFrame(hcat([a.u for a in abm.aout]...)', extract_colnames(abm.aout[1].u)))
        
        mout = DataFrame(hcat([[x.t] for x in abm.mout]...)', [:t]) |> 
        x-> hcat(x, DataFrame(hcat([m.u for m in abm.mout]...)', extract_colnames(abm.mout[1].u)))

        return mout, aout
    end
end

"""
DEBBase Agent.
Each agent owns a parameter collection, state variables and derivatives. 
The states and derivatives include referencs to the global states and derivatives.
"""
mutable struct BaseAgent <: AbstractAgent
    odefuncs::Vector{Function}
    rulefuncs::Vector{Function}
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

        a.odefuncs = Function[
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
        ]

        a.rulefuncs = Function[
            reproduce_opportunistic!,
            die!
        ]

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
