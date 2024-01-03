using Revise
@time using DEBBase
@time using DoseResponse
using Plots, StatsPlots, Plots.Measures
default(leg = false, titlefontsize = 12, legendtitlefontsize = 10, lw = 1.5)
using DataFrames

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
        p = BaseParams(glb = glb, deb = deb)
        isolate_pmoas!(p.deb, [pmoa])
        out_zj = simulator(p)
        out_zj[!,:C_W] .= C_W
        out_zj[!,:pmoa] .= pmoa
        append!(out, out_zj)
    end
end


let pmoas = unique(out.pmoa), num_pmoas = length(pmoas)
    plt = plot(
        layout = (2, num_pmoas), 
        title = hcat(vcat(pmoas, repeat([""], num_pmoas))...), 
        size = (1000,450),
        xlabel = hcat(vcat(
            repeat([""], num_pmoas),
            repeat(["t"], num_pmoas))...), 
        ylabel = ["S" "" "" "" "R" "" "" ""],
        bottommargin = 5mm, leftmargin = 5mm
    )

    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :S, group = :C_W, subplot = j) for (j,pmoa) in enumerate(pmoas)]
    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :R, group = :C_W, subplot = j+num_pmoas) for (j,pmoa) in enumerate(pmoas)]

    display(plt)
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
            drc_params_G = [(1., 4.)],
            drc_params_M = [(1., 4.)],
            drc_params_A = [(1., 4.)],
            drc_params_R = [(1., 4.)],
            drc_params_h = [(1., 4.)]
            )
        p = BaseParams(glb = glb, deb = deb)
        isolate_pmoas!(p.deb, [pmoa])
        out_zj = simulator(p)
        out_zj[!,:C_W] .= C_W
        out_zj[!,:pmoa] .= pmoa
        append!(out, out_zj)
    end
end


let pmoas = unique(out.pmoa), num_pmoas = length(pmoas)
    plt = plot(
        layout = (2, num_pmoas), 
        title = hcat(vcat(pmoas, repeat([""], num_pmoas))...), 
        size = (1000,450),
        xlabel = hcat(vcat(
            repeat([""], num_pmoas),
            repeat(["t"], num_pmoas))...), 
        ylabel = ["S" "" "" "" "R" "" "" ""],
        bottommargin = 5mm, leftmargin = 5mm
    )

    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :S, group = :C_W, subplot = j) for (j,pmoa) in enumerate(pmoas)]
    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :R, group = :C_W, subplot = j+num_pmoas) for (j,pmoa) in enumerate(pmoas)]

    display(plt)
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
            drc_params_G = [(1., 6.)],
            drc_params_M = [(1., 6.)],
            drc_params_A = [(1., 6.)],
            drc_params_R = [(1., 6.)],
            drc_params_h = [(1., 6.)]
            )
        p = BaseParams(glb = glb, deb = deb)
        isolate_pmoas!(p.deb, [pmoa])
        out_zj = simulator(p)
        out_zj[!,:C_W] .= C_W
        out_zj[!,:pmoa] .= pmoa
        append!(out, out_zj)
    end
end


let pmoas = unique(out.pmoa), num_pmoas = length(pmoas)
    plt = plot(
        layout = (2, num_pmoas), 
        title = hcat(vcat(pmoas, repeat([""], num_pmoas))...), 
        size = (1000,450),
        xlabel = hcat(vcat(
            repeat([""], num_pmoas),
            repeat(["t"], num_pmoas))...), 
        ylabel = ["S" "" "" "" "R" "" "" ""],
        bottommargin = 5mm, leftmargin = 5mm
    )

    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :S, group = :C_W, subplot = j) for (j,pmoa) in enumerate(pmoas)]
    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :R, group = :C_W, subplot = j+num_pmoas) for (j,pmoa) in enumerate(pmoas)]

    display(plt)
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
            drc_params_G = [(1., 8.)],
            drc_params_M = [(1., 8.)],
            drc_params_A = [(1., 8.)],
            drc_params_R = [(1., 8.)],
            drc_params_h = [(1., 8.)]
            )
        p = BaseParams(glb = glb, deb = deb)
        isolate_pmoas!(p.deb, [pmoa])
        out_zj = simulator(p)
        out_zj[!,:C_W] .= C_W
        out_zj[!,:pmoa] .= pmoa
        append!(out, out_zj)
    end
end


let pmoas = unique(out.pmoa), num_pmoas = length(pmoas)
    plt = plot(
        layout = (2, num_pmoas), 
        title = hcat(vcat(pmoas, repeat([""], num_pmoas))...), 
        size = (1000,450),
        xlabel = hcat(vcat(
            repeat([""], num_pmoas),
            repeat(["t"], num_pmoas))...), 
        ylabel = ["S" "" "" "" "R" "" "" ""],
        bottommargin = 5mm, leftmargin = 5mm
    )

    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :S, group = :C_W, subplot = j) for (j,pmoa) in enumerate(pmoas)]
    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :R, group = :C_W, subplot = j+num_pmoas) for (j,pmoa) in enumerate(pmoas)]

    display(plt)
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
            drc_params_G = [(1., 10.)],
            drc_params_M = [(1., 10.)],
            drc_params_A = [(1., 10.)],
            drc_params_R = [(1., 10.)],
            drc_params_h = [(1., 10.)]
            )
        p = BaseParams(glb = glb, deb = deb)
        isolate_pmoas!(p.deb, [pmoa])
        out_zj = simulator(p)
        out_zj[!,:C_W] .= C_W
        out_zj[!,:pmoa] .= pmoa
        append!(out, out_zj)
    end
end


let pmoas = unique(out.pmoa), num_pmoas = length(pmoas)
    plt = plot(
        layout = (2, num_pmoas), 
        title = hcat(vcat(pmoas, repeat([""], num_pmoas))...), 
        size = (1000,450),
        xlabel = hcat(vcat(
            repeat([""], num_pmoas),
            repeat(["t"], num_pmoas))...), 
        ylabel = ["S" "" "" "" "R" "" "" ""],
        bottommargin = 5mm, leftmargin = 5mm
    )

    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :S, group = :C_W, subplot = j) for (j,pmoa) in enumerate(pmoas)]
    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :R, group = :C_W, subplot = j+num_pmoas) for (j,pmoa) in enumerate(pmoas)]

    display(plt)
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
            drc_params_G = [(1., 100.)],
            drc_params_M = [(1., 100.)],
            drc_params_A = [(1., 100.)],
            drc_params_R = [(1., 100.)],
            drc_params_h = [(1., 100.)]
            )
        p = BaseParams(glb = glb, deb = deb)
        isolate_pmoas!(p.deb, [pmoa])
        out_zj = simulator(p)
        out_zj[!,:C_W] .= C_W
        out_zj[!,:pmoa] .= pmoa
        append!(out, out_zj)
    end
end


let pmoas = unique(out.pmoa), num_pmoas = length(pmoas)
    plt = plot(
        layout = (2, num_pmoas), 
        title = hcat(vcat(pmoas, repeat([""], num_pmoas))...), 
        size = (1000,450),
        xlabel = hcat(vcat(
            repeat([""], num_pmoas),
            repeat(["t"], num_pmoas))...), 
        ylabel = ["S" "" "" "" "R" "" "" ""],
        bottommargin = 5mm, leftmargin = 5mm
    )

    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :S, group = :C_W, subplot = j) for (j,pmoa) in enumerate(pmoas)]
    [@df out[out.pmoa .== pmoa,:] plot!(plt, :t, :R, group = :C_W, subplot = j+num_pmoas) for (j,pmoa) in enumerate(pmoas)]

    display(plt)
end
