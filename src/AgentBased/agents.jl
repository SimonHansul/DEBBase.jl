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
    global_statevar_indices::Vector{Int}
    
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

        a.global_statevar_indices = findall(x -> x in keys(global_statevars), keys(a.u))

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