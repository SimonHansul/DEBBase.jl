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
    glb = GlobalParams()
    spc_ode = SpeciesParams(
        Z = Dirac(1.), #Truncated(Normal(1, 0.1), 0, Inf),
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
using DataFrames, DataFramesMeta
using Chain

@time sim = AgentBased.simulator(p, saveat = 1) |> agent_record_to_df;
sort!(sim, :t);

sim = combine(groupby(sim, :id)) do df
    df[!,:dI] = diffvec(df.I) ./ diffvec(df.t)
    return df
end;
sim.t = round.(sim.t, sigdigits = 2)


eval_df = @chain sim, simODE begin
    [@select(x, :t, :S, :R, :H) for x in _]
    leftjoin(_[1], _[2], on = :t, makeunique = true)
    drop_na
    @transform(
        :relerr_S = :S ./ :S_1,
        :relerr_H = :H ./ :H_1,
        :relerr_R = :R ./ :R_1, 

        :abserr_S = abs.(:S .- :S_1),
        :abserr_H = abs.(:H .- :H_1),
        :abserr_R = abs.(:R .- :R_1), 
    )
end

@testset begin
    @info("Checking deviation between ODE and ABM simulations")
    @test sum(eval_df.abserr_S .> 15) == 0
    @test sum(eval_df.abserr_R .> 15) == 0
    @test sum(eval_df.abserr_H .> 15) == 0
end


begin
    @df sim plot(
        plot(:t, [:age, :t], color = [1 :gray], ylabel = "age"),
        plot(:t, :X_p, group = :id, ylabel = "X_p"), 
        plot(:t, :f_X, group = :id, ylabel = "f_X", ylim = (0, 1.01)),
        plot(:t, :X_emb, group = :id, ylabel = "X_emb"),
        plot(:t, :dI, group = :id, ylabel = "dI"),
        plot(:t, :S, group = :id, ylabel = "S"),
        xlabel = "t",
        size = (800,500), 
        lw = 2
    )

    @df simODE plot!(:t, :S, subplot = 6)

    hline!([DEBODE.calc_S_max(p.spc)^(2/3) * p.spc.Idot_max_rel_0], subplot = 5, c = :black)
end

using BenchmarkTools
@benchmark AgentBased.simulator(p)
@benchmark DEBODE.simulator(DEBParamCollection())