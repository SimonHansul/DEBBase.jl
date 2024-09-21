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
using DEBBase.Utils

begin # parameter settings
    p = Params()
    p.glb.t_max = 56
    p.glb.N0 = 10
    p.glb.Xdot_in = 30_000
    p.glb.k_V = 0.1 
    p.glb.V_patch = 0.5
    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)


    # setting the TKTD parameters for an initial simulation
    p.spc.k_D_A = [0.5]
    p.spc.e_A = [1.]
    p.spc.b_A = [2.]
    isolate_pmoas!(p.spc, ["A"])
end

@time sim = exposure(
    x -> replicates(AgentBased.simulator, p, 3),
    p,
    [0., 0.5, 1.]
)

using DataFrames
popdata = combine(groupby(sim, [:t, :C_W])) do df
    DataFrame(
        N_tot = nrow(df),
        W_tot = sum(df.S) + sum(df.R)
    )
end

popdata.C_W = [x.C_W[1] for x in eachrow(popdata)]

using Plots, StatsPlots, DEBBase.Figures
@df popdata groupedlineplot(:t, :N_tot, :C_W)
