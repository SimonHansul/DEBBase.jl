
# initialize.jl
# functions to initial states, agents, models, etc.

const X_EMB_INT_REL::Float64 = 0.001 # the (assumed) initial amount of structure, relative to the mass of an egg

"""
    init_substates_agent(p::AbstractParamCollection)
    
Initialize the agent substates, i.e. the state variables of a specific agent.
"""
function init_substates_agent(p::AbstractParamCollection)
    return ComponentArray( # initial states
        X_emb = Float64(p.agn.X_emb_int), # initial mass of vitellus
        S = Float64(p.agn.X_emb_int * X_EMB_INT_REL), # initial structure is a small fraction of initial reserve // mass of vitellus
        H = Float64(0), # maturity
        H_b = Float64(0), # maturity at birth (will be derived from model output)
        R = Float64(0), # reproduction buffer
        I_emb = Float64(0), # ingestion from vitellus
        I_p = Float64(0), # ingestion from external food resource
        I = Float64(0), # total ingestion
        A = Float64(0), # assimilation
        M = Float64(0), # somatic maintenance
        J = Float64(0), # maturity maintenance 

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
    init_substates_global(p::AbstractParamCollection)
Initialize the global substates, i.e. the global state variables such as simulation time.
"""
function init_substates_global(p::AbstractParamCollection)
    return ComponentArray(
            X_p = Float64(p.glb.Xdot_in), # initial resource abundance equal to influx rate
            C_W = (p.glb.C_W) # external stressor concentrations
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
    initialize_statevars!(agent::DEBAgent, abm::DEBABM)::Nothing
Initialize agent-level state variables.
"""
function initialize_statevars!(agent::DEBAgent, abm::DEBABM)::Nothing

    agent.u = ComponentArray(
        glb = Ref(abm.u),
        agn = init_substates_agent(agent.p)
    )
    
    return nothing
end

"""
    initialize_statevars!(abm::DEBABM)
Initialize ABM-level state variables.
"""
function initialize_statevars!(abm::DEBABM)
    abm.u = ComponentArray(
        X_p = Float64(p.glb.Xdot_in), # initial resource abundance equal to influx rate
        C_W = (p.glb.C_W), # external stressor concentrations
    )
end

"""
    init_agents!(abm::DEBABM)::nothing
Initialize the population of agents.
"""
function initialize_agents!(abm::DEBABM)::Nothing

    abm.agents = [] # initialize a vector of agents with undefined values and defined length

    for _ in 1:abm.p.glb.N0 # for the number of initial agents
        push!(abm.agents, abm.p.glb.AgentType(abm)) # initialize an agent and add it to the vector of agents
    end

    return nothing
end