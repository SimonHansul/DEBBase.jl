"""
Set parameters with common prefix to the value of a reference parameter. 
E.g. 

```
set_equal!(spc, :Idot_max_rel, :lrv)
```

sets all parameters whose names start with `Idot_max_rel` equal to `Idot_max_rel_lrv`.
"""
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

Put values from `source` into `sink`. 

- `source::AbstractParams`: A param struct from which we want to extract all parameters.
- `sink::AbstractParams`: All values from `sink` will be assigned to corresponding fields of `sink`.
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
    into!(sink::AbstractParams, source::AbstractParams)

Put values from `source` into `sink`. 

- `source::AbstractParams`: A param struct from which we want to extract all parameters.
- `Sink::DataType`: All values from `sink` will be assigned to a new instance of type `Sink`, assuming that a constructor `Sink()` exists.
"""
function into(Sink::DataType, source::AbstractParams)

    sink_instance = Sink()
    into!(sink_instance, source)

    return sink_instance
end

"""
    trim!(spc::AbstractParams)

Trim TKTD parameters to the smallest number of stressors indicated for any PMoA.
For example `k_D_G = [1., 0.], k_D_M = [1.]` will be trimmed to `k_D_G = [1.], k_D_M = [1.]`.
"""
function trim!(spc::AbstractParams)

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

function isolate_pmoas!(spc::AbstractParams, pmoas::Vector{String}, z::Int64)::Nothing
    spc = isolate_pmoas(spc, pmoas, z)
    return nothing
end

function isolate_pmoas!(spc::AbstractParams, pmoas::Vector{String})::Nothing
    spc = isolate_pmoas(spc, pmoas)
    return nothing
end

"""
    getstat(data::Vector{DataFrame}, statistic::String)

Retrieve the mean of a statistic from a data frame with columns 
"statistic" (listing different kinds of summary statistics) and "value".
"""
function getstat(data::DataFrame, statistic::String)
    return data[data.statistic .== statistic,:].value |> mean
end

function getstat(data::Vector{DataFrame}, statistic::String)
    return getstat(data[2], statistic)
end
