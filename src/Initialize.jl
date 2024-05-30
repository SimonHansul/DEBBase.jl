
# Initialize.jl
# functions to initial states, agents, models, etc. 
# note that many initializations also occur through constructors defined in Structures.jl. here we deal mainly with ComponentArrays

const X_EMB_INT_REL::Float64 = 0.001 # the (assumed) initial amount of structure, relative to the mass of an egg

"""
    init_substates_agent(p::AbstractParamCollection)
    
Initialize the agent substates, i.e. the state variables of a specific agent.
"""
function init_substates_agent(p::AbstractParamCollection)
    return ComponentArray( # initial states
        
    )
end

"""
    init_substates_global(p::AbstractParamCollection)
Initialize the global substates, i.e. the global state variables such as simulation time.
"""
function init_substates_global(p::AbstractParamCollection)
    return ComponentArray(
            X_p = Float64(p.glb.Xdot_in), # initial resource abundance equal to influx rate
            C_W = p.glb.C_W, # external stressor concentrations
            N_tot = p.glb.N0
        )
end

"""
    initialize_statevars(p::AbstractParamCollection, pagnt::ComponentVector{Float64})::ComponentArray 
For initialization of ODE simulator, initialize the component vector of state variables, `u`, based on common parameters `p`.
"""
function initialize_statevars(p::AbstractParamCollection)::ComponentArray 
    return ComponentArray(
        glb = init_substates_global(p),
        agn = init_substates_agent(p)
    )
end

"""
    initialize_statevars!(agent::AbstractAgent)
Initialize agent-level state variables.
"""
function initialize_statevars!(agent::AbstractAgent, abm::AbstractABM)::Nothing
    agent.u = ComponentArray(
        glb = abm.u,
        agn = init_substates_agent(agent.p)
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
        N_tot = p.glb.N0
    )
end

"""
    init_agents!(abm::AbstractABM)::nothing
Initialize the population of 
"""
function initialize_agents!(abm::AbstractABM)::Nothing

    abm.agents = [] # initialize a vector of agents with undefined values and defined length

    for i in 1:abm.p.glb.N0 # for the number of initial agents
        push!(abm.agents, DEBAgent(abm)) # initialize an agent and add it to the vector of agents
    end

    return nothing
end