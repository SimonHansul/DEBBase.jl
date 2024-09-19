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
using DEBBase.AgentBased

# first we extend the species parameters to include additional params needed for the ABM
begin
    glb = GlobalParams(N0 = 1)
    spc_ode = SpeciesParams(
        Z = Truncated(Normal(1, 0.1), 0, Inf),
        K_X_0 = 100 / 0.05
        ) |> ntfromstruct
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
end;

a = DEBAgent(p, DEBODE.initialize_global_statevars(p), 1)
@test a isa AbstractDEBAgent


m = AgentBased.ABM(p);
begin 
    @info "Induce variability in agent params"
    Idot_int = map(x -> x.p.agn.Idot_max_rel_0, m.agents) 
    extremevals = extrema(Idot_int) 
    relrange = extremevals[2]/extremevals[1]
    @test relrange > 0.01
end

@info "Induce variability in initial agent states"
@test map(x -> x.u.S, m.agents) |> extrema |> x -> (x[2]/x[1]) > 0.01

@testset begin
    @info "Transfer global statevars from model to agent"
    @info "Independency of agent statevars from global statevars"
    m.u.C_W = [rand()]
    AgentBased.get_global_statevars!(m.agents[1], m)
    @test m.agents[1].u.C_W == m.u.C_W
    m.agents[1].u.C_W .+= 1
    @test m.agents[1].u.C_W != m.u.C_W
end

@testset begin 
    @info "Update agent statevars using Euler!()"
    m = AgentBased.ABM(p)
    a = m.agents[1]
    S0 = a.u.S
    DEBODE.DEBODE_IA!(a.du, a.u, a.p, m.t)
    AgentBased.Euler!(a.u, a.du, m.dt, Symbol[keys(DEBODE.initialize_agent_statevars(p))...])
    @test a.u.S > S0
end

@test begin
    @info "Run model_step! without err"
    AgentBased.model_step!(m)
    true
end

using Plots
default(leg = false)

@info "Execute model step"
@test begin 
    m = AgentBased.ABM(p)
    AgentBased.model_step!(m)
    true
end

begin
    glb = GlobalParams(
        t_max = 30,
        N0 = 1, 
        Xdot_in = 100, 
        k_V = 0.
        )
    spc_ode = SpeciesParams(
        Z = Truncated(Normal(1, 0.1), 0, Inf),
        K_X_0 = 100 / 0.05,
        Idot_max_rel_0 = 0.
        ) |> ntfromstruct
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
end;

using StatsPlots
using DEBBase.Utils
default(leg = false)

using DEBBase.AgentBased


# FIXME: everything happens too fast...
# dI exceeds maximum

# dX = X_in - k_V * X
# dX -= dI
# X = X - dX


p.glb.Xdot_in = 10.
p.glb.t_max = 3.
dt = 1/24
sim = AgentBased.simulator(p; dt = dt) |> agent_record_to_df

begin
    @df sim plot(
        plot(eachindex(:t), [:age, :t], color = [1 :gray], ylabel = "age"),
        plot(eachindex(:t), :X_p, ylabel = "X_p"), 
        plot(eachindex(:t), :f_X, ylabel = "f_X", ylim = (0, 1.01)),
        plot(eachindex(:t), :X_emb, ylabel = "X_emb"),
        plot(eachindex(:t), diffvec(:I), ylabel = "dI"),
        plot(eachindex(:t), :S, ylabel = "S"),
        xlabel = "step",
        marker = true, size = (800,500)
        )

    hline!([
        DEBODE.calc_S_max(p.spc)^(2/3) * p.spc.Idot_max_rel_0], 
        c = :gray, subplot = 5, 
        leg = true, label = "dI_max_abs"
        )
end
