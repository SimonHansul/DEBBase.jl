using Pkg; Pkg.activate("test")
using Revise

using DEBBase.Utils
using DEBBase.DEBODE
using DEBBase.ABC
using DEBBase.Figures

using Plots, StatsPlots, StatsBase
default(legendtitlefontsize = 10)
using Setfield

# we first load a config file

using YAML
using DataFrames, DataFramesMeta

config = YAML.load_file("test/config/data_config_example.yml")

root = "./"

#=
## Preprocessing of the raw data files
=#

"""
Regression equation for drymass vs length in D. magna.
Source: https://doi.org/10.3390/su14159216 
"""
drymass_from_length(L_mm) = 0.008 * exp(L_mm) - 0.009

begin
    # mapping categorical food treatment label to value in mg dry mass/d
    map_food_levels = Dict(
        "D" => 0.4, 
        "F" => 0.1
    )
    growth = CSV.read("test/data/dp3_test_data_growth_controls_dead.csv", DataFrame)
    growth[!,:drymass] = @. drymass_from_length(growth.Length)
    growth[!,:food_treatment] = [Utils.which_in(x, ["_F_", "_D_"])[2:2] for x in growth.Label]
    growth[!,:food_level] = [map_food_levels[x] for x in growth.food_treatment]

    rename!(
        growth, 
        :t_day => :t_birth,
        :Length => :length,
        )
    select!(growth, :food_level, :t_birth, :length, :replicate, :drymass)

    CSV.write("test/data/dp3_test_data_growth_tidied.csv", growth)

    growth_agg = combine(groupby(growth, [:t_birth, :food_level])) do df
        DataFrame(
            length_mean = mean(skipmissing(df.length)),
            drymass_mean = mean(skipmissing(df.drymass))
        )
    end

    @df growth groupedlineplot(
        :t_birth, :drymass, :food_level, 
        marker = true, mswidth = 2, lw = 2, label = hcat(unique(fround.(:food_level))...), legendtitle = "Food level (mg/d)",
        xlabel = "Time since birth (d)", ylabel = "Dry mass (mg)"
        )
    @df growth scatter!(:t_birth, :drymass, group = :food_level, palette = palette(:default)[1:2], malpha = .75, msize = 3, mswidth = 0.1, label = "")

    @df growth_agg plot(
        :t_birth, :drymass_mean, group = :food_level, 
        c = :black, lw = 2, marker = [:diamond :circle], linestyle = [:dash :solid],
        xlabel = "Time since birth (d)", ylabel = "Mean dry mass (mg)", legendtitle = "Food level (mg/d)"
        )

    CSV.write("test/data/dp3_test_data_growth_tidied_agg.csv", growth_agg)

    # growth stats
    max_drymass(df) = maximum(df.drymass)
    drymass_at_birth(df) = @subset(df, :t_birth .== 0).drymass |> x -> length(x) > 0 ? x[1] : missing

    growth_stats = combine(groupby(growth, [:food_level, :replicate])) do df
        DataFrame(
            max_drymass = max_drymass(df),
            drymass_at_birth = drymass_at_birth(df)
        )
    end

    @df growth_stats boxplot(string.(:food_level), :max_drymass)
end

YAML.load_file("test/data/dpulex_stats.yml")["dpulex_stats"]



using Parameters, DataStructures, DataFrames, CSV

@with_kw struct Dataset
    time_resolved::OrderedDict{String,DataFrame} = OrderedDict()
    scalar::OrderedDict{String,OrderedDict} = OrderedDict()
end

data = Dataset()

for ts_data in config["time_resolved"]
    df = CSV.read(joinpath(root, ts_data["path"]), DataFrame)
    data.time_resolved[ts_data["name"]] = df
end


for sc_data in config["scalar"]
    
    data.scalar = OrderedDict(zip(df.statistic, df.value))
end
data.time_resolved


