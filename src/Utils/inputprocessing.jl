function set_equal!(spc::AbstractParams, prefix::SS, ref_suffix::SS) where SS <: Union{Symbol,String}
    prefix = String(prefix)
    ref_suffix = String(ref_suffix)

    ref_paramname = prefix*"_"*ref_suffix
    ref_param = getproperty(spc, Symbol(ref_paramname))
    paramnames = String[String.(fieldnames(typeof(spc)))...]
    filter!(x -> occursin(String(prefix), x), paramnames)
    filter!(x -> x != String(ref_paramname), paramnames)

    for paramname in paramnames
        setproperty!(spc, Symbol(paramname), ref_param)
    end
end


"""
    into!(sink::AbstractParams, source::AbstractParams)::Nothing

Put values from `source` into `sink`, modifying `sink`.

- `source`: A param struct from which we want to extract all parameters.
- `sink`: All values from `source` will be assigned to corresponding fields of `sink`.
"""
function into!(sink::AbstractParams, source::AbstractParams)::Nothing
    skippedfields = Symbol[]

    for field in fieldnames(typeof(source))
        if field in fieldnames(typeof(sink))
            setproperty!(sink, field, getproperty(source, field))
        else
            push!(skippedfields, field)
        end
    end

    @warn("The following are not fields of $(typeof(sink)) and were skipped: \n $skippedfields")

    return nothing
end


""" 
    into(Sink::DataType, source::AbstractParams)::Sink

Put values from `source` into a new instance of `Sink`.

- `source`: A param struct from which we want to extract all parameters.
- `Sink`: All values from `sink` will be assigned to a new instance of type `Sink`, assuming that a constructor `Sink()` exists.
"""
function into(Sink::DataType, source::AbstractParams)::Sink

    sink_instance = Sink()
    into!(sink_instance, source)

    return sink_instance
end

"""
    trim!(spc::AbstractParams)

Trim TKTD parameters to the smallest number of stressors indicated for any PMoA.
For example `k_D_G = [1., 0.], k_D_M = [1.]` will be trimmed to `k_D_G = [1.], k_D_M = [1.]`.

- `spc`: An instance of a species params type
"""
function trim!(spc::AbstractSpeciesParams)

    num_stressors = min(
        length(spc.k_D_G),
        length(spc.k_D_M),
        length(spc.k_D_A),
        length(spc.k_D_R),
        length(spc.k_D_h)
    )

    spc.k_D_G = spc.k_D_G[1:num_stressors]
    spc.k_D_M = spc.k_D_M[1:num_stressors]
    spc.k_D_A = spc.k_D_A[1:num_stressors]
    spc.k_D_R = spc.k_D_R[1:num_stressors]
    spc.k_D_h = spc.k_D_h[1:num_stressors]

    spc.drc_functs_G = spc.drc_functs_G[1:num_stressors]
    spc.drc_functs_M = spc.drc_functs_M[1:num_stressors]
    spc.drc_functs_A = spc.drc_functs_A[1:num_stressors]
    spc.drc_functs_R = spc.drc_functs_R[1:num_stressors]
    spc.drc_functs_h = spc.drc_functs_h[1:num_stressors]
    
    spc.drc_params_G = spc.drc_params_G[1:num_stressors]
    spc.drc_params_M = spc.drc_params_M[1:num_stressors]
    spc.drc_params_A = spc.drc_params_A[1:num_stressors]
    spc.drc_params_R = spc.drc_params_R[1:num_stressors]
    spc.drc_params_h = spc.drc_params_h[1:num_stressors]
end


"""
    isolate_pmoas(spc::AbstractParams, pmoas::Vector{String}, z::Int64)::AbstractParams

Isolate the indicated PMoAs for chemical stressor `z`. 
That means, turn off all PMoAs (including lethal effects `h`) except for those indicated in `pmoas`-Vector. 
This is done through the toxicokinetic rate constant.

- `spc`: species params
- `pmoas`: Vector of PMoAs we want to keep active
- `z`: Index of chemical stressor we want to be affected

---
## Example


```Julia

p = Params()
p.k_d_G = [1., 1.]
p.k_d_M = [1., 1.]
p.k_d_M = [1., 1.]

p.spc = isolate_pmoas(p, ["G"], 1) # -> for the first stressor, all k_D-values but k_d_G will be 0, second stressor remains untouched

```


"""
function isolate_pmoas(spc::AbstractParams, pmoas::Vector{String}, z::Int64)::AbstractParams
    deactivate = filter(x -> !(x in pmoas), ["G", "M", "A", "R", "h"])
    for j in deactivate
        let fieldname = Symbol("k_D_$(j)")
            k_D = getfield(spc, fieldname)
            k_D[z] = 0.
            setfield!(spc, fieldname, k_D)
        end
    end
    return spc
end

"""
    isolate_pmoas(spc::AbstractParams, pmoas::Vector{String})::AbstractParams

Isolates the given PMoAs for all stressors.
"""
function isolate_pmoas(spc::AbstractParams, pmoas::Vector{String})::AbstractParams
    deactivate = filter(x -> !(x in pmoas), ["G", "M", "A", "R", "h"])
    for j in deactivate
        let fieldname = Symbol("k_D_$(j)")
            k_D = getfield(spc, fieldname)
            k_D .= 0.
            setfield!(spc, fieldname, k_D)
        end
    end
    return spc
end

"""
Mutating version if isolate_pmoas.
"""
function isolate_pmoas!(spc::AbstractParams, pmoas::Vector{String}, z::Int64)::Nothing
    spc = isolate_pmoas(spc, pmoas, z)
    return nothing
end

function isolate_pmoas!(spc::AbstractParams, pmoas::Vector{String})::Nothing
    spc = isolate_pmoas(spc, pmoas)
    return nothing
end

isolate_pmoas!(p::Union{NamedTuple,AbstractParamCollection}, pmoas::Vector{String}) = isolate_pmoas!(p.spc, pmoas)
isolate_pmoas!(p::Union{NamedTuple,AbstractParamCollection}, pmoas::Vector{String}, z::Int64) = isolate_pmoas!(p.spc, pmoas, z)

function getstat(data::DataFrame, statistic::String)
    return data[data.statistic .== statistic,:].value |> mean
end

function getstat(data::Vector{DataFrame}, statistic::String)
    return getstat(data[2], statistic)
end


C2K(T_C::Float64) = T_C + 273.15


abstract type AbstractDataset end

"""
    struct Dataset <: AbstractDataset

Datatype to store calibration data, consisting of `time_resolved` and `scalar_data`. 
This is most useful if data from different sources and/or of different types is pulled 
together for calibration.

What goes into a `Dataset` instance is best defined through a config file 
(cf `test/config/data_config_example.yml`.)

Time-resolved data is assumed to be stored as `csv` file and will be parsed as data frame. 
Scalar data can be stored as `csv` (parsed as data frame) or `yml` (parsed as dict).

Fields:

- time_resolved: Dictionary of time-resolved (tabular) data
- scalar: Dictionary of scalar data (tabular or in dict-form)
- time_vars: Dictionary indicating the time-column for each time-resolved dataset
- grouping_vars: Dictionary indicating additional grouping variables for time-resolved and scalar datasets
- response_vars: Dictionary indicating response variables for time-resolved and scalar datasets
"""
@with_kw struct Dataset <: AbstractDataset
    time_resolved::OrderedDict{String,DataFrame} = OrderedDict()
    scalar::OrderedDict{String,Union{DataFrame,OrderedDict}} = OrderedDict()
    time_vars::Dict = Dict()
    grouping_vars::Dict = Dict("time_resolved" => Dict(), "scalar" => Dict())
    response_vars::Dict = Dict("time_resolved" => Dict(), "scalar" => Dict())
    weights::Dict = Dict("time_resolved" => Dict(), "scalar" => Dict())
end


function read_from_path(path::AbstractString)

    file_extension = split(path, ".")[end]
    @assert file_extension in ["csv", "yml"] "Unknown file extension $(file_extension), expecting csv or yml"

    if file_extension == "csv"
        return CSV.read(path, DataFrame)
    end

    if file_extension == "yml"
        return OrderedDict(YAML.load_file(path))
    end
end

function load_time_resolved_data!(data::AbstractDataset, config::Dict)::Nothing

    for ts_data in config["time_resolved"]
        df = CSV.read(ts_data["path"], DataFrame)
        data.time_resolved[ts_data["name"]] = df
        
        # add entries for response and independent vars

        @assert "grouping_vars" in keys(ts_data) "Independent variables entry missing for time-resolved data in config file"
        @assert "response_vars" in keys(ts_data) "Response variables entry missing for time-resolved data in config file"
        
        data.time_vars[ts_data["name"]] = Symbol(ts_data["time_var"])

        data.grouping_vars["time_resolved"][ts_data["name"]] = Symbol.(
            ts_data["grouping_vars"]
        )

        data.response_vars["time_resolved"][ts_data["name"]] = Symbol.(
            ts_data["response_vars"]
        )

        data.weights["time_resolved"][ts_data["name"]] = ts_data["weight"]
    end

    return nothing
end

function load_scalar_data!(data::AbstractDataset, config::Dict)::Nothing
    for sc_data in config["scalar"]

        data.scalar[sc_data["name"]] = read_from_path(sc_data["path"])
        
        # grouping_vars are optional here
        @assert "response_vars" in keys(sc_data) "Response variables entry missing for scalar data in config file"

        if "grouping_vars" in keys(sc_data)
            data.grouping_vars["scalar"][sc_data["name"]] = Symbol.(
                sc_data["grouping_vars"]
            )
        else
            data.grouping_vars["scalar"][sc_data["name"]] = Symbol[]
        end

        data.response_vars["scalar"][sc_data["name"]] = Symbol.(
            sc_data["response_vars"]
        )

        data.weights["scalar"][sc_data["name"]] = sc_data["weight"]
    end

    return nothing
end


"""
    data_from_config(path_to_config::AbstractString; datatype::DataType = Dataset)::datatype

Load a dataset from data config file. 

```Julia
usng DEBBase.Utils
data = Utils.data_from_config("path/to/config.yml")
```
"""
function data_from_config(path_to_config::AbstractString; datatype::DataType = Dataset)::datatype

    config = YAML.load_file(path_to_config)
    data = datatype()
    load_time_resolved_data!(data, config)
    load_scalar_data!(data, config)

    return data
end

