begin
    using Pkg; Pkg.activate("tests")
    using Revise
    @time using DEBBase
    @time using DoseResponse
    using Plots, StatsPlots, Plots.Measures
    default(leg = false, titlefontsize = 12, legendtitlefontsize = 10, lw = 1.5)
    using DataFrames, DataFramesMeta
    using StatsBase
    const TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"") 
end

#=
Simulate single stressors with different PMoAs
=#
begin   
    out = DataFrame()
    for pmoa in ["G", "M", "A", "R"]
        for C_W in round.(10 .^ range(log10(0.1), log10(1.), length = 5), sigdigits = 2)
            glb = GlobalBaseParams(t_max = 42., C_W = [C_W])
            deb = DEBBaseParams(
                k_D_G = [10.], 
                k_D_M = [10.], 
                k_D_A = [10.], 
                k_D_R = [10.], 
                k_D_h = [10.], 
                drc_params_G = [(1., 2.)],
                drc_params_M = [(1., 2.)],
                drc_params_A = [(1., 2.)],
                drc_params_R = [(1., 2.)],
                drc_params_h = [(1., 2.)]
                )
            p = BaseParamCollection(glb = glb, deb = deb)
            isolate_pmoas!(p.deb, [pmoa])
            out_zj = DEBBase.simulator(p)
            out_zj[!,:C_W] .= C_W
            out_zj[!,:pmoa] .= pmoa
            append!(out, out_zj)
        end
    end
    out = DEBBase.relative_response(out, [:S, :R], :C_W; groupby_vars = [:t, :pmoa])
end

out = DataFrame()
for pmoa in ["G", "M", "A", "R"]
    for C_W in round.(10 .^ range(log10(0.1), log10(2.), length = 5), sigdigits = 2)
        glb = GlobalBaseParams(t_max = 42., C_W = [C_W])
        deb = DEBBaseParams(
            k_D_G = [1.], 
            k_D_M = [1.], 
            k_D_A = [1.], 
            k_D_R = [1.], 
            k_D_h = [1.], 
            drc_params_G = [(1., 2.)],
            drc_params_M = [(1., 2.)],
            drc_params_A = [(1., 2.)],
            drc_params_R = [(1., 2.)],
            drc_params_h = [(1., 2.)]
            )
        p = BaseParamCollection(glb = glb, deb = deb)
        isolate_pmoas!(p.deb, [pmoa])
        out_zj = DEBBase.simulator(p)
        out_zj[!,:C_W] .= C_W
        out_zj[!,:pmoa] .= pmoa
        append!(out, out_zj)
    end
end
