abstract type AbstractDEBAgent end

CAUSE_OF_DEATH = Dict(
    0 => "none",
    1 => "age"
)

@with_kw mutable struct DEBAgent <: AbstractDEBAgent
    id::Int
    age::Real
    cause_of_death::Int
    du::ComponentVector
    u::ComponentVector
    p::Union{AbstractParamCollection,NamedTuple}

    time_since_last_repro::Real
    cum_offspring::Int64
    cohort::Int64
    
    function DEBAgent(p::Union{AbstractParamCollection,NamedTuple}, global_statevars::ComponentVector, id; cohort = 0)
        a = new() # create empty agent instance

        a.p = ( # agent holds its own parameter object
            glb = p.glb, # global params
            spc = p.spc, # species params
            agn = (; # agent params
                ntfromstruct(ODEAgentParams(p.spc))...,
                (
                    a_max = rand(p.spc.a_max)
                )
            ) # agn
        ) # a.p
        
        # vcat does not seem to work on more than two component arrays (returns Vector instead), hence the pipe syntax
        a.u = vcat(
            global_statevars,
            initialize_agent_statevars(a.p)
        )

        a.du = similar(a.u)
        a.du .= 0.

        a.id = id
        a.cohort = cohort
        a.age = 0.
        a.cause_of_death = 0
        a.time_since_last_repro = 0.
        a.cum_offspring = 0.
        
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
    saveat::Real
    agent_record::Vector{ComponentVector}
    global_statevar_names::Vector{Symbol}
    global_statevar_indices::Vector{Int64}
    agent_statevar_names::Vector{Symbol}
    agent_statevar_indices::Vector{Int64}

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
    function ABM(p::Union{AbstractParamCollection,NamedTuple}; dt = 1/24, saveat = 1)::ABM
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
        m.agent_record = ComponentVector[]

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

get_recorded_agent_var_indices(m::AbstractDEBABM) = map(x -> x in m.recorded_agent_vars,  keys(a.u)) |> BitVector

function get_global_statevars!(a::AbstractDEBAgent, m::AbstractDEBABM)::Nothing
    # TODO: replacte setproperty!/getproperty with direct indexing
    for var in keys(m.u)
        setproperty!(a.u, var, getproperty(m.u, var))
        setproperty!(a.du, var, getproperty(m.du, var))
    end
end

function set_global_statevars!(m::AbstractDEBABM, a::AbstractDEBAgent)::Nothing
    for var in keys(m.u)
        setproperty!(m.u, var, getproperty(a.u, var))
        setproperty!(m.du, var, getproperty(a.du, var))
    end
end

function agent_step_rulebased!(a::AbstractDEBAgent, m::AbstractDEBABM)::Nothing
    a.age += m.dt 

    if condition_juvenile(a.u, m.t, a) <= 0
        effect_juvenile!(a)
    end

    if condition_adult(a.u, m.t, a) <= 0
        effect_adult!(a)
    end

    if a.age >= a.p.agn.a_max
        a.cause_of_death = 1
    end
    
    if (a.u.S / a.u.S_max_hist) < 0.5
        if rand() < exp(-0.7 * m.dt)
            a.cause_of_death = 2
        end
    end

    #if (a.u.S / a.u.S_max_hist) < 0.5 # TODO:make the threshold a parameter
    #    if rand() < exp(-0.7 * m.dt) 
    #        a.u.cause_of_death = 2
    #    end
    #end

    if a.time_since_last_repro >= a.p.spc.tau_R
        let num_offspring = trunc(a.u.R / a.u.X_emb_int)
            for _ in 1:num_offspring
                m.idcount += 1
                push!(m.agents, DEBAgent(a.p, m.u, m.idcount; cohort = a.cohort + 1))
                a.u.R -= a.u.X_emb_int
            end
            a.time_since_last_repro = 0.
        end
    else
        a.time_since_last_repro += m.dt
    end

    return nothing
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

    DEBODE_agent_IA!(a.du, a.u, a.p, m.t)

    Euler!(a.u, a.du, m.dt, m.agent_statevar_indices)
    agent_step_rulebased!(a, m)
    set_global_statevars!(m, a)

    return nothing
end

function Euler!(u::ComponentVector, du::ComponentVector, dt::Real, statevar_indices::Vector{Int})::Nothing
    
    u[statevar_indices] .+= du[statevar_indices] .* dt
    
    #for ui in statevar_names
    #    setproperty!(u, ui, getproperty(u, ui) + getproperty(du, ui) * dt)
    #end

    return nothing
end

function record_agent!(a::AbstractDEBAgent, m::AbstractDEBABM)::Nothing

    if isapprox(m.t % m.saveat, 0, atol = m.dt)
    push!(
        m.agent_record,
        vcat(
            ComponentVector(
                t = m.t, 
                id = a.id, 
                age = a.age, 
                cause_of_death = a.cause_of_death
                ),
            a.u
        ) # vcat
    ) # push
    end

    return nothing
end

filter_agents!(m) = m.agents = filter(x -> x.cause_of_death == 0, m.agents)

"""
    step_all_agents!(m::AbstractDEBABM)::Nothing

Executes agent_step! for all agents in the model. 
Records agent states.
Filters agents vector to remove dead indiviualds. 
Agents which die will be recorded a last time before they are removed.
"""
function step_all_agents!(m::AbstractDEBABM)::Nothing
    for a in m.agents
        agent_step!(a, m)
        record_agent!(a, m)
    end
    filter_agents!(m)

    return nothing
end

function model_step!(m::AbstractDEBABM)::Nothing
    # calculate global derivatives
    # change in resource abundance, chemical stressor exposure etc.
    
    DEBODE_global!(m.du, m.u, m.p, m.t)
    step_all_agents!(m)
    
    # global statevars are updated after agent derivatives are calculated
    # this is important because agents affect global states using mutating operators
    Euler!(m.u, m.du, m.dt, m.global_statevar_indices) 
    m.u.X_p = max(0, m.u.X_p) # HOTFIX : negative resource abundances cause chaos

    m.t += m.dt

    return nothing
end

function simulator(
    p::Union{NamedTuple,AbstractParamCollection}; 
    dt = 1/24, 
    saveat = 1
    )
    #@info "Running ABM simulation with t_max=$(p.glb.t_max)"
    
    m = ABM(p; dt = dt, saveat = saveat)

    while !(m.t > m.p.glb.t_max)
        #isapprox(m.t % 14, 0, atol = m.dt) ? @info("t=$(m.t)") : nothing

        model_step!(m)
    end

    return m
end

function agent_record_to_df(
    m::AbstractDEBABM; 
    #cols::Vector{Symbol} = [:t, :id, :embryo, :juvenile, :adult, :age, :f_X, :S, :I, :A, :M, :H, :R, :X_p, :X_emb]
    )::DataFrame
    cols = vcat(
        [:t, :id],
        [keys(m.agent_record[1])...]
    )  |> unique
    
    hcat([map(x -> getproperty(x, y), m.agent_record) for y in cols]...) |> 
    x -> DataFrame(x, cols)
end
