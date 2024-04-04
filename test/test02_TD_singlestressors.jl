if occursin("terminal", abspath(PROGRAM_FILE))
    @info("Loading packages")
    using Pkg; Pkg.activate("test")
    using Plots, StatsPlots, Plots.Measures
    using SHUtils
    using DataFrames, DataFramesMeta
    using ProgressMeter
    default(leg = false, lw = 1.5)
    TAG = splitpath(@__FILE__)[end] |> x -> split(x, ".")[1] |> String    
    using DEBFigures

    using Revise
    using DEBBase
end
TAG = replace(splitpath(@__FILE__)[end], ".jl" =>"")


#=
Simulate single stressors with different PMoAs
=#
begin   
    out = DataFrame()
    for pmoa in ["G", "M", "A", "R"]
        for C_W in round.(10 .^ range(log10(0.1), log10(1.), length = 5), sigdigits = 2)
            glb = GlobalParams(t_max = 42., C_W = [C_W])
            spc = SpeciesParams(
                k_D_G = [10.], 
                k_D_M = [10.], 
                k_D_A = [10.], 
                k_D_R = [10.], 
                k_D_h = [10.], 
                e_G = [1.],
                e_M = [1.],
                e_A = [1.],
                e_R = [1.],
                e_h = [1.],
                b_G = [1.],
                b_M = [1.],
                b_A = [1.],
                b_R = [1.],
                b_h = [1.],
                )
            p = DEBParamCollection(glb = glb, spc = spc)
            isolate_pmoas!(p.spc, [pmoa])
            out_zj = simulator(p)
            out_zj[!,:C_W] .= C_W
            out_zj[!,:pmoa] .= pmoa
            append!(out, out_zj)
        end
    end
    out = DEBBase.relative_response(out, [:S, :R], :C_W; groupby_vars = [:t, :pmoa])

    @df out plot(
        groupedlineplot(:t, :y_S, :C_W)
    )
end

