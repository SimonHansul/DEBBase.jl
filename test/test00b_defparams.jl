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

    using SHUtils
    using Revise
    @time using DEBBase    
end

#=


=#
@testset begin # testing the default parameters
    theta = DEBParamCollection()
    theta.glb.t_max = 56.
    theta.spc.Z = Dirac(1.)
    yhat = DEBBase.simulator(theta)
    @df yhat plot(
        plot(:t, :S),
        plot(:t, :H)
     ) |> display

    @test isapprox(maximum(yhat.H), theta.spc.H_p, rtol = 1e-2)
    @test isapprox(maximum(yhat.S), DEBBase.calc_S_max(theta.spc), rtol = 0.1)
end;

 
#=
Basic test of @replicates macro

=#

@testset begin
    theta = DEBParamCollection()
    theta.spc.Z = Truncated(Normal(1., 0.1), 0, Inf)
    yhat = @replicates DEBBase.simulator(theta) 10

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

    @test cvs.S > 0.05
    @test cvs.H > 0.05
    @test cvs.R > 0.05
end;


#=
Basic test of @compose macro. <br>
Here we simulate again the default parameters, but use the @compose macro to define the ODE system.
=#

#@testset begin
#
#    functions = [ # define ODE system as list of derivative functions
#        DEBBase.Idot!,
#        DEBBase.Adot!,
#        DEBBase.Mdot!,
#        DEBBase.Jdot!,
#        DEBBase.Sdot!,
#        DEBBase.Hdot!,
#        DEBBase.H_bdot!,
#        DEBBase.Rdot!,
#        DEBBase.X_pdot!,
#        DEBBase.X_embdot!,
#        DEBBase.Ddot!,
#        DEBBase.C_Wdot!
#    ]
#
#    system! = DEBBase.@compose functions # use @compose to put the functions together
#    
#    theta = DEBParamCollection()
#    theta.glb.t_max = 56.
#
#    @time yhat = DEBBase.simulator(theta; system = system!)
#
#    plt = @df yhat plot(
#        plot(:t, :S, ylabel = "S"),
#        plot(:t, :H, ylabel = "H"), 
#        plot(:t, :R, ylabel = "R"),
#        plot(:t, diffvec(:I) ./ diffvec(:t), ylabel = "Idot"), 
#        title = "@compose ", titlefontsize = 10,
#        xlabel = "t"
#    )
#
#    display(plt)
#    
#    @test isapprox(maximum(yhat.H), theta.spc.H_p, rtol = 1e-2)
#    @test isapprox(maximum(yhat.S), DEBBase.calc_S_max(theta.spc), rtol = 0.1)
#end

theta = DEBParamCollection()

using DEBFigures

theta = DEBParamCollection()
theta.spc.Z = Truncated(Normal(1, 0.05), 0, Inf)

@time yhat = sweep(
    :(@replicates simulator(theta) 10),
    theta.spc, :Idot_max_rel, 
    range(10, 25, length = 10)
    )

@df yhat groupedlineplot(
    :t, :S, :Idot_max_rel, 
    xlabel = "Time (d)", ylabel = "Structure (μgC)", leftmargin = 5mm,
    palette = palette([:teal, :purple], 10), 
    leg = :outertopright, 
    label = hcat(fround.(unique(:Idot_max_rel))...), legtitle = "Idot_max_rel"
    )



