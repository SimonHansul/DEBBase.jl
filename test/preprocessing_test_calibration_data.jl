# preprocessing of the data used for testing calibration routines
# this script is not used for the tests themselves

using Pkg; Pkg.activate("test")

using Plots, StatsPlots, StatsBase
default(legendtitlefontsize = 10)
using DataFrames, DataFramesMeta, CSV, YAML

using Revise

using DEBBase.Utils
using DEBBase.DEBODE
using DEBBase.ABC
using DEBBase.Figures

#=
## Preprocessing of the raw data files
=#

"""
Regression equation for drymass vs length in D. magna.
Source: https://doi.org/10.3390/su14159216 
"""
drymass_from_length(L_mm) = 0.008 * exp(L_mm) - 0.009

using CSV
begin # tidying the growth data
    # mapping categorical food treatment label to value in mg dry mass/d
    map_food_levels = Dict(
        "_Co_" => 0.4, 
        "_F_" => 0.1
    )
    growth = CSV.read("test/data/dp3_test_data_growth_controls_dead.csv", DataFrame)
    growth[!,:drymass] = @. drymass_from_length(growth.Length)
    growth[!,:food_treatment] = [Utils.which_in(x, ["_Co_", "_F_"]) for x in growth.Label]
    growth[!,:food_level] = [map_food_levels[x] for x in growth.food_treatment]

    rename!(
        growth, 
        :t_day => :t_birth,
        :Length => :length,
        )
    select!(growth, :food_level, :t_birth, :length, :replicate, :drymass)

    CSV.write("test/data/dp3_test_data_growth_tidy.csv", growth)
end

begin # aggregating the tidy growth data
    growth_agg = combine(groupby(growth, [:t_birth, :food_level])) do df
        DataFrame(
            length_mean = mean(skipmissing(df.length)),
            length_sd = std(skipmissing(df.length)),
            length_n = sum(@. !ismissing(df.length)),
            drymass_mean = mean(skipmissing(df.drymass)),
            drymass_sd = std(skipmissing(df.drymass)),
            drymass_n = sum(@. !ismissing(df.drymass))
        )
    end
    CSV.write("test/data/dp3_test_data_growth_tidy_agg.csv", growth_agg)

    @df growth groupedlineplot(
        :t_birth, :drymass, :food_level, 
        marker = true, mswidth = 2, lw = 2, label = hcat(unique(fround.(:food_level))...), legendtitle = "Food level (mg/d)",
        xlabel = "Time since birth (d)", ylabel = "Dry mass (mg)"
        )
    @df growth scatter!(:t_birth, :drymass, group = :food_level, palette = palette(:default)[1:2], malpha = .75, msize = 3, mswidth = 0.1, label = "")

    @df growth_agg plot(
        :t_birth, :drymass_mean, group = :food_level, 
        yerr = :drymass_sd,
        c = :black, lw = 2, marker = [:diamond :circle], linestyle = [:dash :solid],
        xlabel = "Time since birth (d)", ylabel = "Mean dry mass (mg)", legendtitle = "Food level (mg/d)"
        )

end

begin # calculating additional summary statistics from the growth data
    t_max = maximum(growth.t_birth)
    # growth stats
    max_drymass(df) = df[df.t_birth .== t_max, :drymass] |> x -> length(x)>0 ? x[1] : missing
    drymass_at_birth(df) = @subset(df, :t_birth .== 0).drymass |> x -> length(x) > 0 ? x[1] : missing

    growth_stats = combine(groupby(growth, [:food_level, :replicate])) do df
        DataFrame(
            max_drymass = max_drymass(df),
            drymass_at_birth = drymass_at_birth(df)
        )
    end

    CSV.write("test/data/dp3_test_data_growth_stats.csv", growth_stats)

    plot(
        (@df growth_stats[:,[:food_level,:max_drymass]] |> drop_na boxplot(string.(:food_level), :max_drymass)),
        (@df growth_stats[:,[:food_level,:drymass_at_birth]] |> drop_na boxplot(string.(:food_level), :drymass_at_birth)),
        leg = false, 
        xlabel = "Food level (mg/d)", ylabel = ["Maximum dry mass (mg)" "Dry mass at birth (mg)"], 
        ylim = [(0.025, 0.12) (0, 0.015)], 
        boxwidth = 2
    )
end

begin # aggregating the growth stat da
    growth_stats_agg = combine(groupby(growth_stats, :food_level)) do df
        DataFrame(
            max_drymass_mean = Utils.robustmean(df.max_drymass),
            max_drymass_sd = std(skipmissing(df.max_drymass)),
            max_drymass_n = sum(@. !ismissing(df.max_drymass)),
            drymass_at_birth_mean = Utils.robustmean(df.drymass_at_birth),
            drymass_at_birth_sd = std(skipmissing(df.drymass_at_birth)),
            drymass_at_birth_n = sum(@. !ismissing(df.drymass_at_birth))
        )
    end
    CSV.write("./test/data/dp3_test_data_growth_stats_agg.csv", growth_stats_agg)
end

begin # tidying the repro data
    map_food_levels = Dict(
        "F" => 0.1,
        "C" => 0.4
    )
    repro = CSV.read("test/data/dp3_test_data_reproduction_controls_dead.csv", DataFrame)
    repro[!,:food_treatment] = [Utils.which_in(x, ["F-", "Co-"])[1:1] for x in repro.trtm]
    repro[!,:food_level] = [map_food_levels[uppercase(x)] for x in repro.food_treatment]

    rename!(
        repro,
        :rep => :replicate,
        :t_day => :t_birth
    )

    select!(repro, :t_birth, :replicate, :food_level, :cum_repro)

    CSV.write("test/data/dp3_test_data_repro_tidy.csv", repro)
    
    @df repro groupedlineplot(
        :t_birth, :cum_repro, :food_level, marker = true, lw = 2, 
        xlabel = "Time since birth (d)", ylabel = "Cumulative reproduction [#]", 
        label = hcat(unique(:food_level)...), legendtitle = "Food level (mg/d)"
    )
end

begin # aggregating the reproduction data
    repro_agg = combine(groupby(repro, [:t_birth, :food_level])) do df
        DataFrame(
            cum_repro_mean = Utils.robustmean(df.cum_repro),
            cum_repro_sd = std(skipmissing(df.cum_repro)),
            cum_repro_n = sum(@. !ismissing(df.cum_repro))
        )
    end

    CSV.write("test/data/dp3_test_data_repro_tidy_agg.csv", repro_agg)

    @df repro_agg plot(
        :t_birth, :cum_repro_mean, group = :food_level, 
        yerr = :cum_repro_sd,
        c = :black, lw = 2, marker = [:diamond :circle], linestyle = [:dash :solid],
        xlabel = "Time since birth (d)", ylabel = "Mean cumulative \n reproduction (#)", 
        legendtitle = "Food level (mg/d)"
        )
end

begin # computing additional summary statistics from repro data
    """
    Time between birth and first hatching of offspring
    """
    birth_to_first_hatch(x) = x[x.cum_repro .> 0,:].t_birth |> x-> length(x)>0 ? x[1] : missing

    """
    Final cumulative reproduction 
    """
    final_cum_repro(x, t_max) = x[x.t_birth.==t_max,:].cum_repro |> x-> length(x)>0 ? x[1] : missing

    repro_stats = combine(groupby(repro, [:replicate, :food_level])) do df
        DataFrame(
            birth_to_first_hatch = birth_to_first_hatch(df),
            final_cum_repro = final_cum_repro(df, t_max)
        )
    end

    CSV.write("test/data/dp3_test_data_repro_stats.csv", repro_stats)

    plot(
        (@df drop_na(repro_stats[:,[:food_level,:birth_to_first_hatch]]) boxplot(string.(:food_level), :birth_to_first_hatch)),
        (@df drop_na(repro_stats[:,[:food_level,:final_cum_repro]]) boxplot(string.(:food_level), :final_cum_repro)), 
        xlabel = "Food level (mg/d)", 
        ylabel = ["Birth to first hatch (d)" "Final cumulative reproduction (#)"], leg = false, 
        ylim = [(4, 10) (20, 100)]
    )
end

begin # aggregating the repro summary statistics
    repro_stats_agg = combine(groupby(repro_stats, [:food_level])) do df
        DataFrame(
            birth_to_first_hatch_mean = Utils.robustmean(df.birth_to_first_hatch),
            birth_to_first_hatch_sd = std(skipmissing(df.birth_to_first_hatch)),
            birth_to_first_hatch_n = sum(@. !ismissing(df.birth_to_first_hatch)),

            final_cum_repro_mean = Utils.robustmean(df.final_cum_repro),
            final_cum_repro_sd = std(skipmissing(df.final_cum_repro)),
            final_cum_repro_n = sum(@. !ismissing(df.final_cum_repro))
        )
    end

    CSV.write("test/data/dp3_test_data_repro_stats_agg.csv", repro_stats_agg)
end
