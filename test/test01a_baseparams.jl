#=
Testing for changes in baseline parameters

We conduct parameter sweeps to test whether growth, maturation and reproduction respond to changes in baseline parameters as expected.
=#

@info("Loading packages")
@time begin
    using Pkg; Pkg.activate("test")
    using Plots, StatsPlots, Plots.Measures
    using SHUtils
    using DataFrames, DataFramesMeta
    using ProgressMeter
    default(leg = false, lw = 1.5)
    TAG = splitpath(@__FILE__)[end] |> x -> split(x, ".")[1] |> String    

    using Revise
    using DEBBase
end

# TODO: implement baseparam test as in AmphiDEB


begin # debugging for small species
    @info("Setting up parameters")
    ref = DEBBaseParams() # reference params
    Smax_ref = DEBBase.calc_S_max(ref)

    deb = copy(ref) # adjusted params
    deb.Idot_max_rel = 0.1 # 
    deb.Idot_max_rel_emb = deb.Idot_max_rel
    deb.eta_AS = 0.9

    Smax_deb = DEBBase.calc_S_max(deb)
    Z = Smax_deb / Smax_ref
    deb.X_emb_int = ref.X_emb_int * Z

    glb = GlobalBaseParams()
    glb.Xdot_in *= Z

    println(
        "initial reserves: $(round(deb.X_emb_int, sigdigits = 3)) \n",
        "expected maximum structure: $(round(Smax_deb, sigdigits = 3))"
        )

    #=
    Why no growth during larval phase for low Idot_max_rel?

    Observations:
    - S stays way below theoretical S_max
    - Idot_p stays constant


    Problem?

    - f(X) -> no
    - Idot_p -> stays constant, makes sense
    - Adot -> stays constant, makes sense
    - Sdot -> goes to 0, does not make sense!!!
    - MCov -> satys constant and >0, makes sense


    So the problem has to be in the calculation of Sdot for larvae?



    =#


    @info("Running simulations")
    y = DEBBase.simulator(
        BaseParamCollection(glb = glb, deb = deb), 
        dt = 1/280
        )

    @df y plot(
        plot(:t, :S, ylabel = "S"),
        plot(:t, (:X_p ./ 0.05) ./ ((:X_p ./ 0.05) .+ deb.K_X), ylim = (0,1.01), ylabel = "f(X)"),
        plot(:t, :X_p, ylabel = "X_p"),
        plot(:t, diffvec(:M) ./ diffvec(:t), ylabel = "Mdot"),
        plot(:t, diffvec(:A) ./ diffvec(:t), ylabel = "Adot"), 
        plot(:t, diffvec(:S) ./ diffvec(:t), ylabel = "Sdot"),
        plot(:t, (deb.kappa .* (diffvec(:A) ./ diffvec(:t)) .- diffvec(:M) ./ diffvec(:t)), ylabel = "Mcov", ylim = (0, Inf)),
        xlabel = "t", 
        size = (800,500), lw = 2, 
        leftmargin = 5mm
    )
    hline!([Smax_deb], color = :gray, linestyle = :dash, subplot = 1)
end


