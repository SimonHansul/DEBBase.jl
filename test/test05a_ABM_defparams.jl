# script to test the DEB-ABM on the default parameter set
using Pkg; Pkg.activate("test")
using Parameters
using NamedTupleTools
using Distributions
using ComponentArrays
using Test

using Revise 
@time using DEBBase.DEBODE
using DEBBase.ParamStructs
using DEBBase.DoseResponse
#using DEBBase.DEBABM

# first we extend the species parameters to include additional params needed for the ABM
begin
    glb = GlobalParams(N0 = 3)
    spc_ode = SpeciesParams(Z = Truncated(Normal(1, 0.1), 0, Inf)) |> ntfromstruct
    agn_ode = DEBODE.ODEAgentParams(DEBParamCollection().spc) |> ntfromstruct

    spc = (; 
        spc_ode...,
        (
            a_max = Truncated(Normal(60, 6), 0, Inf), # maximum life span
            tau_R = 2. # reproduction interval
        )...
    )

    agn = (;
        agn_ode...,
        (
            a_max = rand(spc.a_max)
        )
    )

    p = (glb = glb, spc = spc, agn = agn)
end

a = DEBAgent(p, DEBODE.initialize_global_statevars(p), 1)
@test a isa AbstractDEBAgent

m = ABM(p)
@info "Induce variability in agent params"
@test map(x -> x.p.agn.Idot_max_rel_0, m.agents) |> extrema |> x-> (x[2]/x[1])>0.01
@info "Induce variability in initial agent states"
@test map(x -> x.u.S, m.agents) |> extrema |> x -> (x[2]/x[1]) > 0.01


@testset begin
    @info "Transfer global statevars from model to agent"
    @info "Independency of agent statevars from global statevars"
    m.u.C_W = [rand()]
    get_global_statevars!(m.agents[1], m)
    @test m.agents[1].u.C_W == m.u.C_W
    m.agents[1].u.C_W .+= 1
    @test m.agents[1].u.C_W != m.u.C_W
end

@testset begin 
    @info "Update agent statevars using Euler!()"
    m = ABM(p)
    a = m.agents[1]
    S0 = a.u.S
    DEBODE.DEBODE_IA!(a.du, a.u, a.p, m.t)
    Euler!(a.u, a.du, m.dt)
    @test a.u.S > S0
end

@test begin
    @info "Run model_step! without err"
    model_step!(m)
    true
end

using Plots
@info "Execute model step"
@test begin 
    m = ABM(p)
    model_step!(m)
    true
end

sim = simulator(p);

p.glb.t_max

t = map(x -> x.t, sim.agent_record)
id = map(x -> x.id, sim.agent_record)
age = map(x -> x.age, sim.agent_record)
S = map(x -> x.S, sim.agent_record)
cod = map(x -> x.cause_of_death, sim.agent_record)

plot(
    plot(t, S, group = id), 
    plot(t, age, group = id),
    leg = false
)

