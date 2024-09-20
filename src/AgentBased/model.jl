abstract type AbstractDEBABM end

mutable struct ABM <: AbstractDEBABM
    agents::Vector{DEBAgent}
    du::ComponentVector
    u::ComponentVector
    p::Union{NamedTuple,AbstractParamCollection}
    t::Real
    dt::Real
    idcount::Int
    saveat::Real
    agent_record::Vector{ComponentVector}
    global_statevar_names::Vector{Symbol}
    global_statevar_indices::Vector{Int64}
    agent_statevar_names::Vector{Symbol}
    agent_statevar_indices::Vector{Int64}

    record_agents::Bool

    """
        DEBABM(p::Union{AbstractParamCollection,NamedTuple}; dt = 1/24)::DEBABM
    
    Initialization of the DEBABM structure.
    
    args
        - `p`: A compatible parameter collection (`ABMParamCollection` or correspondingly structured named tuple)
    
    kwargs
        - `dt`: Model time step [t]
    """
    function ABM(p::Union{AbstractParamCollection,NamedTuple}; dt = 1/24, saveat = 1, record_agents::Bool)::ABM
        m = new()
        m.agents = Vector{DEBAgent}(undef, p.glb.N0)
        m.u = initialize_global_statevars(p)
        m.global_statevar_names = Symbol[keys(m.u)...]
        m.global_statevar_indices = collect(eachindex(m.global_statevar_names))
        m.du = similar(m.u) 
        m.p = p
        m.t = 0
        m.dt = dt
        m.saveat = saveat
        m.idcount = 0

        m.record_agents = record_agents
        m.agent_record = ComponentVector[]

        # initialize agents
        for i in 1:p.glb.N0
            m.idcount += 1
            m.agents[i] = DEBAgent(p, m.u, m.idcount)
        end

        m.agent_statevar_names = Symbol[keys(initialize_agent_statevars(p))...]

        if p.glb.N0 > 0
            m.agent_statevar_indices = findall(x -> x in m.agent_statevar_names, keys(m.agents[1].u))
        else
            m.agent_statevar_indices = Int64[]
        end

        return m
    end
end