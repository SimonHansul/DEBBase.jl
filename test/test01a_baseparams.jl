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

#=
debugging for small species (S_max ≈ 1e-5)...
we set the issues with very small species aside for a second and focus on the size range we realistically need. 
expected egg mass for our species is in the range of 1 mg. 
applying some margin of error, we can set X_emb_int to 0.1 mg and use this to zoom Idot_max_rel, resulting in the smallest species we may need to simulate.
This is a species with maximum Size around 1 mg.
=#
begin 
    @info("Setting up parameters")
    ref = DEBBaseParams() # reference params
    Smax_ref = DEBBase.calc_S_max(ref) # S_max of the reference species

    deb = copy(ref) # adjusted params
    deb.X_emb_int = 1e-4 # egg mass is fixed 
    Z = deb.X_emb_int / ref.X_emb_int # zoom factor based on egg sizes

    deb.Idot_max_rel *= Z^(1/3) # zoomed ingestion rate
    deb.Idot_max_rel_emb = deb.Idot_max_rel # embryoinc ingestion rate equal to non-embryonic
    deb.K_X *= Z^(1/3) # zoomed half-saturation constant (K_X ∝ Idot_max_rel)

    Smax_deb = DEBBase.calc_S_max(deb) # maximum size of the reference species

    glb = GlobalBaseParams()
    glb.Xdot_in *= Z
    glb.t_max = 30.

    println(
        "zoom factor: $(round(Z, sigdigits = 2)) \n",
        "initial reserves: $(round(deb.X_emb_int, sigdigits = 4)) \n",
        "expected maximum structure: $(round(Smax_deb, sigdigits = 3))"
        )
    

    @info("Running simulations")
    y = DEBBase.simulator(
        BaseParamCollection(glb = glb, deb = deb)
        )

    @df y plot(
        plot(:t, :S, ylabel = "S"),
        plot(:t, (:X_p ./ 0.05) ./ ((:X_p ./ 0.05) .+ deb.K_X), ylim = (0,1.01), ylabel = "f(X)"),
        plot(:t, :X_p, ylabel = "X_p"),
        plot(:t, diffvec(:I) ./ diffvec(:t), ylabel = "Idot"),
        plot(:t, diffvec(:A) ./ diffvec(:t), ylabel = "Adot"), 
        plot(:t, diffvec(:M) ./ diffvec(:t), ylabel = "Mdot"),
        plot(:t, diffvec(:S) ./ diffvec(:t), ylabel = "Sdot"),
        plot(:t, (deb.kappa .* (diffvec(:A) ./ diffvec(:t)) .- diffvec(:M) ./ diffvec(:t)), ylabel = "Mcov", ylim = (0, Inf)),
        xlabel = "t", 
        size = (800,500), lw = 2, 
        leftmargin = 5mm
    )

    hline!([Smax_deb], color = :gray, linestyle = :dash, subplot = 1)
    hline!([deb.Idot_max_rel * Smax_deb^(2/3)], color = :gray, linestyle = :dash, subplot = 4)
end
