# script to test the DEB-ABM on the default parameter set 
# the most critical test here is a comparison between the ODE and ABM output

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
using DEBBase.DEBABM

# first we extend the species parameters to include additional params needed for the ABM
begin
    glb = GlobalParams(N0 = 1)
    spc_ode = SpeciesParams(
        Z = Truncated(Normal(1, 0.1), 0, Inf),
        K_X_0 = 100 / 0.05
        ) |> ntfromstruct
    agn_ode = DEBODE.AgentParams(Params().spc) |> ntfromstruct

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

a = DEBAgent(p, DEBODE.initialize_global_statevars(p), 1, AgentParams)
@test a isa AbstractDEBAgent
m = DEBABM.ABM(p);
@test m isa AbstractDEBABM

@testset begin 
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
    DEBABM.get_global_statevars!(m.agents[1], m)
    @test m.agents[1].u.C_W == m.u.C_W
    m.agents[1].u.C_W .+= 1
    @test m.agents[1].u.C_W != m.u.C_W
end

@testset begin 
    @info "Update agent statevars using Euler!()"
    m = DEBABM.ABM(p)
    a = m.agents[1]
    S0 = a.u.S
    DEBODE.DEBODE_IA!(a.du, a.u, a.p, m.t)
    DEBABM.Euler!(a.u, a.du, m.dt, collect(eachindex(a.u)))
    @test a.u.S > S0
end

@test begin
    @info "Run model_step! without err"
    DEBABM.model_step!(m)
    true
end

using Plots
default(leg = false)

@info "Execute model step"
@test begin 
    m = DEBABM.ABM(p)
    DEBABM.model_step!(m)
    true
end

begin # setting up parameters
    glb = GlobalParams(
        Xdot_in = 2400
    )
    spc_ode = SpeciesParams(
        Z = Dirac(1.), #Truncated(Normal(1, 0.1), 0, Inf),
        ) |> ntfromstruct
    agn_ode = DEBODE.AgentParams(Params().spc) |> ntfromstruct

    spc = (; 
        spc_ode...,
        (

            a_max = Truncated(Normal(60, 6), 0, Inf), # maximum life span
            tau_R = Inf # reproduction interval
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

using DEBBase.DEBABM
using DataFrames, DataFramesMeta
using Chain

sim = DEBABM.simulator(p) 
simODE = DEBODE.simulator(Params(spc = SpeciesParams(tau_R = Inf)))

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
    # We accept an absolute error of 0.1 μgC for an individual that can reach >300 μgC in structure
    @test sum(eval_df.abserr_S .> 0.1) == 0
    @test sum(eval_df.abserr_R .> 0.1) == 0
    @test sum(eval_df.abserr_H .> 0.1) == 0
end

begin
    plt = @df sim plot(
        plot(:t, :S, group = :id, ylabel = "S", linestyle = :dash),
        plot(:t, :H, group = :id, ylabel = "H", linestyle = :dash, leg = true, label = "ABM"),
        plot(:t, :R, group = :id, ylabel = "R", linestyle = :dash),
        xlabel = "t",
        size = (800,500), 
        lw = 2
    )

    @df simODE plot!(:t, :S, subplot = 1)
    @df simODE plot!(:t, :H, subplot = 2, label = "ODE")
    @df simODE plot!(:t, :R, subplot = 3)

    #hline!([DEBODE.calc_S_max(p.spc)^(2/3) * p.spc.Idot_max_rel_0], subplot = 5)
end

using BenchmarkTools
@info "Comparing simulation times of ODE and ABM: We expect ODE to be faster than ABM"
b = @benchmark DEBABM.simulator(p)
bmedian_agentbased = median(b.times)
bODE = @benchmark DEBODE.simulator(Params())
bmedian_ODE = median(bODE.times)
@info("ODE is $(round(bmedian_agentbased/bmedian_ODE, sigdigits = 2)) times faster than ABM") 

using DEBBase.DEBABM
begin
    glb = GlobalParams(
        t_max = 365,
        N0 = 10,
        Xdot_in = 30_000,
        k_V = 0.1,
        V_patch = 0.5
        )
    spc_ode = SpeciesParams(
        Z = Truncated(Normal(1, 0.1), 0, Inf),
        K_X_0 = 10_000
        ) |> ntfromstruct
    agn_ode = DEBODE.AgentParams(Params().spc) |> ntfromstruct

    spc = (; 
        spc_ode...,
        (
            a_max = Truncated(Normal(60, 6), 0, Inf), # maximum life span
            tau_R = 2.5 # reproduction interval
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

using DEBBase.Figures
using DEBBase.DEBODE

p.glb.t_max = 56.

begin
    @info "Generating output plot for ABM with defaults"
    @time "Computation time using 10 threaded replicates of the ABM with t_max=$(p.glb.t_max):" sim = treplicates(
        DEBABM.simulator, p, 3; saveat = 1
        )
    sort!(sim, :t);

    sim_agg = combine(groupby(sim, [:t, :replicate, :cohort])) do df
        DataFrame(
            N_tot = nrow(df),
            W_tot = sum(df.S) + sum(df.R),
            S_mean = mean(df.S),
            R_mean = mean(df.R),
            fX_mean = mean(df.f_X),
            Xp_mean = mean(df.X_p)
        )
    end

    using StatsBase
    sim.cause_of_death |> countmap

    plot( # show a trajectory for each cohort
        [
            groupedlineplot(sim_agg.t, sim_agg[:,y], sim_agg.cohort, ylabel = y) 
            for y in [:N_tot, :W_tot, :S_mean, :R_mean, :fX_mean, :Xp_mean]]...,
                size = (800,500),
                xlabel = "t", 
        )
end


