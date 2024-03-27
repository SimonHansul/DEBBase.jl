#=
# Testing the basic functionality of macros
=#

begin
    using Pkg; Pkg.activate("test")

    using Plots, StatsPlots
    using Distributions
    using DataFrames
    using Test
    default(leg = false)
    using DifferentialEquations

    using SHUtils
    using Revise
    @time using DEBBase    
end

#=
Basic test of @replicates macro
=#
@test begin
    defparams = BaseParamCollection()
    defparams.deb.Z = Dirac(1.)
    yhat = DEBBase.simulator(defparams)
    @df yhat plot(:t, :S)

    theta = BaseParamCollection()
    theta.deb.Z = Truncated(Normal(1., 0.1), 0, Inf)
    yhat = @replicates DEBBase.simulator(theta) 10

    @df yhat plot(:t, :S, group = :replicate, color = 1)
    true # test is passed if code executes without error
end
#=
Basic test of @compose macro
=#
@test begin
    modelfunctions = [
        DEBBase.Idot!,
        DEBBase.Adot!,
        DEBBase.Mdot!,
        DEBBase.Jdot!,
        DEBBase.Sdot!,
        DEBBase.Hdot!,
        DEBBase.H_bdot!,
        DEBBase.Rdot!,
        DEBBase.X_pdot!,
        DEBBase.X_embdot!,
        DEBBase.Ddot!,
        DEBBase.C_Wdot!
    ]

    testmodel! = @compose modelfunctions
    
    function testsim(pcmn::Ref{BaseParamCollection}; 
        alg = Tsit5(),
        saveat = 1, 
        reltol = 1e-6,
        kwargs...)

        DEBBase.assert!(pcmn)
        pown = DEBBase.initialize_pown()
        DEBBase.agent_variability!(pown, pcmn)
        u = DEBBase.initialize_statevars(pcmn, pown)
        prob = ODEProblem(testmodel!, u, (0, pcmn.x.glb.t_max), (pcmn, pown)) # define the problem
        sol = solve(prob, alg; saveat = saveat, reltol = reltol, kwargs...) # get solution to the IVP
        simout = DEBBase.sol_to_df(sol) # convert solution to dataframe
    
        return simout
    end

    theta = BaseParamCollection()
    theta.glb.t_max = 56.
    tref = Ref(theta)
    global yhat = testsim(tref)

    plt = @df yhat plot(
        plot(:t, :S, ylabel = "S"),
        plot(:t, :H, ylabel = "H"), 
        plot(:t, :R, ylabel = "R"),
        plot(:t, diffvec(:I) ./ diffvec(:t), ylabel = "Idot"), 
        title = "@compose ", titlefontsize = 10,
        xlabel = "t"
    )

    display(plt)
    
    c1 = isapprox(yhat.S[end], DEBBase.calc_S_max(theta.deb), rtol = 0.1)
    c2 = isapprox(yhat.H[end], theta.deb.H_p, rtol = 1e-3)

    c1 & c2
end
