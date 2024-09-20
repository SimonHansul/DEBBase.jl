get_recorded_agent_var_indices(m::AbstractDEBABM) = map(x -> x in m.recorded_agent_vars,  keys(a.u)) |> BitVector


"""
    get_global_statevars!(a::AbstractDEBAgent, m::AbstractDEBABM)::Nothing


"""
function get_global_statevars!(a::AbstractDEBAgent, m::AbstractDEBABM)::Nothing
    
    a.du[a.global_statevar_indices] .= m.du
    a.u[a.global_statevar_indices] .= m.u
    
    return nothing
end


"""
    set_global_statevars!(m::AbstractDEBABM, a::AbstractDEBAgent)::Nothing

Updates global state variables and derivatives according to modifications by an agent.
"""
function set_global_statevars!(m::AbstractDEBABM, a::AbstractDEBAgent)::Nothing

    m.du .= a.du[a.global_statevar_indices]
    m.u .= a.u[a.global_statevar_indices]

    return nothing
end


"""
    agent_step_rulebased!(a::AbstractDEBAgent, m::AbstractDEBABM)::Nothing

Defines the rule-based portion of an agent step. 

The life stage callbacks defined in the DEBODE are used to apply rules for life stage transitions.

A crude rule for starvation mortality is implemented, applying a constant hazard rqte 
when they lose a given fraction of their structural mass. 

Reproduction is assumed to occur in fixed time intervals, according to `spc.tau_R`.
"""
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
    
    let s_f = sig(
        a.u.f_X, a.p.spc.f_Xthr, 
        (1 - a.p.spc.s_min) * a.u.f_X / a.p.spc.f_Xthr + a.p.spc.s_min,
        1.)^m.dt

        if rand() > s_f
            a.cause_of_death = 2.
        end
    end

    #if (a.u.S / a.u.S_max_hist) < 0.66
    #    if rand() < exp(-0.7 * m.dt)
    #        a.cause_of_death = 2
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
At the minimum, this will include life stage transitions, reproduction and death of agents. 
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


"""
    Euler!(u::ComponentVector, du::ComponentVector, dt::Real, statevar_indices::Vector{Int})::Nothing

Applying the Euler scheme to state variables at given indices.

The `statevar_indices` are needed because agents might hold a (reference to) global state variables, 
which we don't wish to modify while updating the agent's state.

args

- `u`: State variables
- `du`: Derivatives
- `dt`: Timestep
- `statevar_indices`: Indices of state variables and derivatives which should be included

"""
function Euler!(u::ComponentVector, du::ComponentVector, dt::Real, statevar_indices::Vector{Int})::Nothing
    
    u[statevar_indices] .+= du[statevar_indices] .* dt

    return nothing
end

"""
    record_agent!(a::AbstractDEBAgent, m::AbstractDEBABM)::Nothing


Pushes the agent state variables to `m.agent_record` in fixed time intervals, according to `m.saveat`.
"""
function record_agent!(a::AbstractDEBAgent, m::AbstractDEBABM)::Nothing

    if m.record_agents && isapprox(m.t % m.saveat, 0, atol = m.dt)
        push!(
            m.agent_record,
            vcat(
                ComponentVector(
                    t = m.t, 
                    id = a.id, 
                    cohort = a.cohort,
                    age = a.age, 
                    cause_of_death = a.cause_of_death
                    ),
                a.u
            ) # vcat
        ) # push
    end

    return nothing
end


"""
    filter_agents!(m::AbstractDEBABM)

Remove dead agents from a model. 
Agents are flagged to be removed if `cause_of_death>0`.
"""
filter_agents!(m::AbstractDEBABM) = m.agents = filter(x -> x.cause_of_death == 0, m.agents)

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