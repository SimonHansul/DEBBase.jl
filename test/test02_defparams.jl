#=
Testing the default parameters

- is the maximum maturity equal to maturity at puberty?
- does the maximum structure match with what is calculated from parameters?
- does the mass balance check out?
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

    @test isapprox(maximum(yhat.H), p.spc.H_p_0, rtol = 1e-2) # test for maximum maturity
    @test isapprox(maximum(yhat.S), DEBODE.calc_S_max(p.spc), rtol = 0.1)

    # calculating cumulative energy budgets to assess mass balance

    
    budget = DataFrame()
    for y in [:S, :M, :J, :R, :Q]
        investment = @subset(yhat, :t .== maximum(:t))[1,y]
        append!(budget, DataFrame(y = y, investment = investment))
    end
    sort!(budget, :y)
    bar(budget.investment, xticks = (1:nrow(budget), budget.y))
    mass_balance = sum(budget.investment) ./ @subset(yhat, :t .== maximum(:t))[1,:t]
    @test isapprox(mass_balance, 1, atol = 0.98)
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
