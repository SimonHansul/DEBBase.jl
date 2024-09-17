#=
Simulate single stressors with different PMoAs
=#


using DEBBase.ParamStructs

"""
    exposure(simcall::Expr, C_Wvec::Vector{Float64}; kwargs...)

Simulate exposure to a single stressor over a Vector of constant exposure concentrations `C_Wvec`. 

"""
function exposure(simulator::Function, params::Union{AbstractParamCollection,NamedTuple}, C_Wvec::Vector{Float64})
    
    let C_W_int = params.glb.C_W # we will modify this value and then reset to the initial value
        sim = DataFrame()

        for C_W in C_Wvec
            params.glb.C_W[1] = C_W
            sim_i = simulator(params)
            append!(sim, sim_i)
        end

        params.glb.C_W = C_W_int 

        return sim
    end
end

using DEBBase.DEBODE

ENV["JULIA_DEBUG"] = Main

params = DEBParamCollection()
isolate_pmoas!(params, ["M"])
params.spc.k_D_M = [1.]
params.spc.e_M = [1.]

sim = exposure(
    simulator, 
    params, 
    [1.]
)

@df sim plot(:t, :D_M_1)

@df sim plot(
    plot(:t, :S, group = :C_W_1, layout = (1,3)),
    plot(:t, :D_M_1, group = :C_W_1, layout = (1,3)), 
    layout = (2,1)
    )


let C_Wvec =  vcat([0], round.(10 .^ range(log10(0.1), log10(1.), length = 5), sigdigits = 2))
    global sim = DataFrame()
    pmoas = ["G", "M", "A", "R"]
    for pmoa in pmoas
        for C_W in C_Wvec
            glb = GlobalParams(t_max = 42., C_W = [C_W])
            spc = SpeciesParams(
                kappa_0 = 0.538,
                k_D_G = [0.5], 
                k_D_M = [0.5], 
                k_D_A = [0.5], 
                k_D_R = [0.5], 
                k_D_h = [0.5], 
                e_G = [1.],
                e_M = [1.],
                e_A = [1.],
                e_R = [1.],
                e_h = [1.],
                b_G = [2.],
                b_M = [2.],
                b_A = [2.],
                b_R = [2.],
                b_h = [2.]
                )
            theta = DEBParamCollection(glb = glb, spc = spc)
            isolate_pmoas!(theta.spc, [pmoa])
            sim_zj = simulator(theta)
            sim_zj[!,:C_W] .= C_W
            sim_zj[!,:pmoa] .= pmoa
            append!(sim, sim_zj)
        end # for C_W in ...
    end # for pmoa in ...
    
    sim = relative_response(sim, [:S, :R], :C_W; groupby_vars = [:t, :pmoa])

    @debug(@subset(sim, :C_W_1 .== 0)[:, [:t, :C_W_1, :R]])
    
    #rankcor = combine(groupby(sim, :pmoa)) do sim_j
    #    @chain sim_j begin
    #        combine(groupby(_, :C_W_1), :y_R => last) 
    #        @aside begin
    #            @debug(@subset(_, ismissing.(:y_R_last)))
    #        end
    #        corspearman(_.C_W_1, _.y_R_last)
    #        return DataFrame(r = _)
    #    end 
    #end

    #@test unique(rankcor.r) == [-1]

    plt = plot(
        layout = (2,4), title = hcat(pmoas...), 
        ylim = (0, 1.01),
        xlabel = "t", ylabel = hcat([gridylabel("y_S", 1, 4), gridylabel("y_R", 1, 4)]...), 
        size = (800,450), 
        bottommargin = 5mm, leftmargin = 2.5mm, lw = 2,
        leg = [false false false true false false false false],
        legtitle = "C_W / e", legendtitlefontsize = 10
        )

    for (j,pmoa) in enumerate(pmoas)
        @df @subset(sim, :pmoa .== pmoa) groupedlineplot!(plt, :t, :y_S, :C_W, subplot = j, label = hcat(unique(:C_W)...))
        @df @subset(sim, :pmoa .== pmoa) groupedlineplot!(plt, :t, :y_R, :C_W, subplot = 4+j, label = hcat(unique(:C_W)...))
    end

    display(plt)
end



@chain sim begin
    @select(:t, :R, :C_W_1, :pmoa, :y_R)
    @subset(:C_W_1 .== 0)
    plot
end

sort!(sim, :t)
@df @subset(sim, :pmoa .== "G") groupedlineplot(:t, :R, :C_W_1)

