using Pkg; Pkg.activate("test")
using Revise
@time using DEBBase.DEBODE, DataFrames, Distributions
using Plots, StatsPlots


p = Params() #
sim = DataFrame()
let Tvec = [17.5, 20., 22.5] # simulate three ambient temperatures
    for T in Tvec
        p.glb.T = 273.15 + T # apply temperature in Kelvin
        sim_i = DEBODE.simulator(p) # generate the predidction
        sim_i[!,:T] .= T
        append!(sim, sim_i) #
    end
end

sort!(sim, :t)

fig1 = @df sim plot(
    plot(:t, :S, group = :T, ylabel = "Structural mass [μg C]", legendtitle = "T [°C]"),
    plot(:t, :R ./ p.spc.X_emb_int_0, group = :T, ylabel = "Cumulative reproduction [#]", leg = false), 
    xlabel = "Age [d]", lw = 1.5
)

savefig(plot(fig1, dpi = 300), "paper/fig1.png")


using DEBBase.AgentBased, Distributions
using Base.Threads
using ProgressMeter
 
# setting the parameters to useful values for population 

p = Params()
p.glb.t_max = 150
p.glb.N0 = 10
p.glb.Xdot_in = 30_000
p.glb.k_V = 0.1 
p.glb.V_patch = 0.5
p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)

Tvec = [17.5, 20., 22.5] # simulate three ambient temperatures
sim = DataFrame()
@showprogress for (i,T) in enumerate(Tvec)  
    p.glb.T = 273.15 + T # apply temperature in Kelvin =
    sim_i = @replicates AgentBased.simulator(p) 10 # generate the predidction
    sim_i[!,:T] .= p.glb.T # add column indicating temperature
    append!(sim, sim_i)
end

popdata = combine(groupby(sim, [:t, :replicate, :T])) do df
    DataFrame(
        N_tot = nrow(df),
        W_tot = sum(df.S) + sum(df.R),
        X_p_median = median(df.X_p),
        S_median = median(df.S),
        freq_emb = sum(df.embryo) ./ nrow(df),
        freq_juv = sum(df.juvenile) ./ nrow(df),
        freq_ad = sum(df.adult) ./ nrow(df),
        ratio_juv = sum(df.juvenile) ./ sum(df.adult)
    )
end

using DataFramesMeta, StatsPlots, Plots.Measures
using DEBBase.Figures

sort!(popdata, :t)

begin
    fig2 = plot(xlabel = "time [d]", legendtitle = "Temperature [°C]", layout = (1,3), size = (1000,500), bottommargin = 5mm, leftmargin = 5mm)

    for (i,T) in enumerate(Tvec .+ 273.15)
        @df @subset(popdata, :T .== T) plot!(
            fig2, subplot = 1,
            :t, :N_tot, group = :replicate, color = i, 
            xlabel = "Time [d]", 
            ylabel = "Total population size [#]", 
            linealpha = 0.5, 
            label = "", title = "Abundance"
            )

        @df @subset(popdata, :T .== T) plot!(
            fig2, subplot = 2,
            :t, :W_tot, group = :replicate, color = i, 
            xlabel = "Time [d]", 
            ylabel = "Total population mass [μg C]", 
            linealpha = 0.5, 
            label = "", title = "Biomass"
            )

        @df @subset(popdata, :T .== T) plot!(
            fig2, subplot = 3,
            :t, :freq_juv, group = :replicate, color = i, 
            xlabel = "Time [d]", 
            ylabel = "Amount of resource [μg C]", 
            linealpha = 0.5, 
            label = "", title = "Resource abundance"
            )
    end

    @df popdata groupedlineplot!(
        fig2, 
        subplot = 1, 
        :t, :N_tot, :T, 
        lw = 2, fillalpha = .25, palette = palette(:default)[1:length(Tvec)], 
        label = hcat(unique(:T) .- 273.15...), leg = true
        )

        
    @df popdata groupedlineplot!(
        fig2, 
        subplot = 2, 
        :t, :W_tot, :T, 
        lw = 2, fillalpha = .25, palette = palette(:default)[1:length(Tvec)], 
        label = hcat(unique(:T) .- 273.15...), leg = false
        )

    @df popdata groupedlineplot!(
        fig2, 
        subplot = 3, 
        :t, :X_p_median, :T,
        lw = 2, fillalpha = .25, palette = palette(:default)[1:length(Tvec)], 
        label = hcat(unique(:T) .- 273.15...), leg = false
        )    


    savefig(plot(fig2, dpi = 300), "paper/fig2.png")
    fig2
end
