#=
# Testing the basic functionality based on default parameters
=#

begin
    using Pkg; Pkg.activate("test")

    using Plots, StatsPlots, Plots.Measures
    default(leg = false)
    using Distributions
    using DataFrames
    using Test
    using OrdinaryDiffEq
    using Chain

    using Revise
    @time using DEBBase.DEBODE
    @time using DEBBase.Utils
end

#=
Testing the default parameters
=#
@testset begin 
    p = DEBParamCollection()
    p.glb.t_max = 56.
    p.spc.Z = Dirac(1.)
    yhat = simulator(p)
    @df yhat plot(
        plot(:t, :S),
        plot(:t, :H)
     ) |> display

    @test isapprox(maximum(yhat.H), p.spc.H_p, rtol = 1e-2) # test for maximum maturity
    @test isapprox(maximum(yhat.S), DEBODE.calc_S_max(p.spc), rtol = 0.1)
end;

#=
Basic test of @replicates macro
=#

@testset begin
    p = DEBParamCollection()
    p.spc.Z = Truncated(Normal(1., 0.1), 0, Inf)
    yhat = @replicates simulator(p) 10

    plt = @df yhat plot(
        plot(:t, :S, group = :replicate, color = 1),
        plot(:t, :H, group = :replicate, color = 1)
    )

    display(plt)

    cvs = @chain yhat begin # compute coefficients of variation in final values
        groupby(:replicate)
        combine(_) do df
            return DataFrame(
                S_max = maximum(df.S),
                H_max = maximum(df.H),
                R_max = maximum(df.R)
            )
        end
        (
            S = std(_.S_max) / mean(_.S_max),
            H = std(_.H_max) / mean(_.H_max),
            R = std(_.R_max) / mean(_.R_max)
        )
    end

    @test cvs.S > 0.05 # test for plausible coefficients of variation in the final values
    @test cvs.H > 0.05
    @test cvs.R > 0.05
end;
