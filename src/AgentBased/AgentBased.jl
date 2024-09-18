
import DEBBase.DEBODE: simulator

abstract type AbstractDEBAgent end

const CAUSE_OF_DEATH = Dict(
    0 => "none",
    1 => "age"
)


@with_kw mutable struct DEBAgent <: AbstractDEBAgent
    id::Int
    du::ComponentVector
    u::ComponentVector
    p::Union{AbstractParamCollection,NamedTuple}
    
    function DEBAgent(p::Union{AbstractParamCollection,NamedTuple}, global_statevars::ComponentVector, id)
        a = new() # create empty agent instance
        a.id = id
        a.p = ( # agent holds its own parameter object
            glb = p.glb, # global params
            spc = p.spc, # species params
            agn = (; # agent params
                ntfromstruct(DEBODE.ODEAgentParams(p.spc))...,
                (
                    a_max = rand(spc.a_max)
                )
            ) # agn
        ) # a.p
        
        # vcat does not seem to work on more than two component arrays (returns Vector instead), hence the pipe syntax
        a.u = vcat(
            global_statevars,
            DEBODE.initialize_agent_statevars(a.p)
        ) |> x-> vcat(x, ComponentVector(age = 0, cause_of_death = 0))

        a.du = similar(a.u)
    
        return a
    end
end


abstract type AbstractDEBABM end

mutable struct ABM <: AbstractDEBABM
    agents::Vector{DEBAgent}
    du::ComponentVector
    u::ComponentVector
    p::Union{NamedTuple,AbstractParamCollection}
    t::Real
    dt::Real
    idcount::Int
    agent_record::Vector{ComponentVector}
    #recorded_agent_vars::Vector{Symbol}
    #recorded_agent_var_indices::BitVector

    """
        DEBABM(p::Union{AbstractParamCollection,NamedTuple}; dt = 1/24)::DEBABM
    
    Initialization of the DEBABM structure.
    
    args
        - `p`: A compatible parameter collection (`ABMParamCollection` or correspondingly structured named tuple)
    
    kwargs
        - `dt`: Model time step [t]
    """
    function ABM(p::Union{AbstractParamCollection,NamedTuple}; dt = 1/24)::ABM
        m = new()
        m.agents = Vector{DEBAgent}(undef, p.glb.N0)
        m.u = DEBODE.initialize_global_statevars(p)
        m.du = similar(m.u) 
        m.p = p
        m.t = 0
        m.dt = dt
        m.idcount = 0
        m.agent_record = ComponentVector[]

        for i in 1:p.glb.N0
            m.idcount += 1
            m.agents[i] = DEBAgent(p, m.u, m.idcount)
        end

        return m
    end
end

get_recorded_agent_var_indices(m::AbstractDEBABM) = map(x -> x in m.recorded_agent_vars,  keys(a.u)) |> BitVector

function get_global_statevars!(a::AbstractDEBAgent, m::AbstractDEBABM)::Nothing
    for var in keys(m.u)
        setproperty!(a.u, var, getproperty(m.u, var))
    end
end

function set_global_statevars!(m::AbstractDEBABM, a::AbstractDEBAgent)::Nothing
    for var in keys(m.u)
        setproperty!(m.u, var, getproperty(a.u, var))
    end
end

function agent_step_rulebased!(a::AbstractDEBAgent, m::AbstractDEBABM)
    a.u.age += m.dt 

    if a.u.age >= a.p.agn.a_max
        a.u.cause_of_death = 1
    end
end


"""
    agent_step!(a::Agent, m::Model)::Nothing

The agent step follows a generic pattern:

First the ODE-portion of the model is executed, and the corresponding state variables are updated using the Euler scheme. 
Then the rule-based portion of the model is executed. These are all the functions which cannot / should not be expressed as part of an ODE.
At the minimum, this will be reproduction and death of agents. 
For a spatially explicit model, movement should also most likely be part of the rule-based portion, 
as well as functions which require direct information exchange between agents.
"""
function agent_step!(a::AbstractDEBAgent, m::AbstractDEBABM)
    get_global_statevars!(a, m)
    DEBODE.DEBODE_IA!(a.du, a.u, a.p, m.t)
    Euler!(a.u, a.du, m.t)
    agent_step_rulebased!(a, m)
    set_global_statevars!(m, a)

    return nothing
end

function Euler!(u::ComponentVector, du::ComponentVector, dt::Real)::Nothing
    for var in keys(u)
        setproperty!(u, var, getproperty(u, var) + getproperty(u, var) * dt)
    end
end

function agent_record!(a::AbstractDEBAgent, m::AbstractDEBABM)::Nothing
    push!(
        m.agent_record,
        vcat(
            ComponentVector(t = m.t, id = a.id),
            a.u
        )
    )

    return nothing
end

filter_agents!(m) = m.agents = filter(x -> x.u.cause_of_death == 0, m.agents)

"""
    agent_step!(m::AbstractDEBABM)

Executes agent_step! for all agents in the model. 
Records agent states.
Filters agents vector to remove dead indiviualds. 
Agents which die will be recorded a last time before they are removed.
"""
function agent_step!(m::AbstractDEBABM)::Nothing
    for a in m.agents
        agent_step!(a, m)
        agent_record!(a, m)
    end
    filter_agents!(m)

    return nothing
end

function model_step!(m::AbstractDEBABM)::Nothing
    DEBODE.DEBODE_global!(m.du, m.u, m.p, m.t)
    Euler!(m.u, m.du, m.t)
    agent_step!(m)

    m.t += m.dt

    return nothing
end


function simulator(p::Union{NamedTuple,AbstractParamCollection}; dt = 1/24)
    m = ABM(p; dt = dt)

    while !(m.t > m.p.glb.t_max)
        model_step!(m)
        println(m.t)
    end

    return m
end