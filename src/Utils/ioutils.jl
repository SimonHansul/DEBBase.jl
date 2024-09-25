function extract_colnames(c::R, k::Symbol) where R <: Real
    return k
end

function extract_colnames(c::AbstractVector, k::Symbol)
    return [Symbol("$(k)_$(i)") for i in 1:length(c)]
end

function extract_colnames(u::ComponentVector)
    colnames = []
    for k in keys(u)
        push!(colnames, extract_colnames(u[k], k))
    end
    return vcat(colnames...)
end

extract_colnames(sol::ODESolution)::Vector{Symbol} = vcat([:t], extract_colnames(sol.u[1]))

ymlparse(x::String) = Meta.parse(x) |> eval
ymlparse(x::Vector{String}) = @. Meta.parse(x) |> eval
ymlparse(x::Vector{N}) where N <: Number = x

"""
    params_from_config(ParamType::DataType, path_to_config::String)::ParamType

Load a configuration file and return an instance of type `P`, which is assumed to be a parameter collection.

- `ParamType`: The type of parameter structure type which is expected to be stored in the config file (e.g. `Params`)
- `path_to_config`: Location of the (YAML) config file`

---

## Example: 

```Julia
stored_params = params_from_config(Params, "config/my_config_file.yml")
```

"""
function params_from_config(ParamType::DataType, path_to_config::String)::ParamType
    config = YAML.load_file(path_to_config)
    p = ParamType()

    # iterate over sub-structures mentioned in config file
    for substruct in keys(config)
        # iterate over parameters mentioned within substructure
        for par in keys(config[substruct])
            let val
                # evaluate the entry, which might or might not be parsed correctly by YAML
                strval = eval(:($config[$substruct][$par]))

                # if the eval returned a numeric type, we don't have to do anything
                if strval isa Number
                    val = strval
                # otherwise, we still have to parse and eval a second time
                else
                    val = ymlparse(strval)
                end

                # some convoluted meta-programming to parse the value and assign it to the param struct...
                :($p.$(Symbol(substruct)).$(Symbol(par)) = $val) |> eval
            end # let val
        end
    end

    return p
end

deserialize(x::Any) = x # most types don't need extra treatment for deserialization

# distributions need extra methods for proper deserialization, 
function deserialize(d::Distribution)
    disttype = "$(typeof(d))"
    distparams = fieldnames(typeof(d))
    deserialized_dist_params = [deserialize(getproperty(d, p)) for p in distparams]
    return "$disttype($(["$(x), " for x in deserialized_dist_params]...))"
end

# truncated distributions are again a special case
function deserialize(d::Truncated)
    untruncated = deserialize(d.untruncated)
    return "Truncated($untruncated, $(d.lower), $(d.upper))"
end


function deserialize_dict!(dict)::Nothing
    for key in keys(dict)
        dict[key] = deserialize(dict[key])
    end
    return nothing 
end


"""
    save_config(path_to_config::String, p::AbstractParamCollection)::Nothing


Save a parameter collection to config file.

- `path_to_config`: location to save the config file (with extension .yml)
- `p`: Instance of an abstract parameter collection holding the parameter values to store

---

## Example

```Julia
using DEBBase.DEBODE, DEBBase.DoseResponse, DEBBase.Utils

# load the default parameters
p = Params() 

# make some modifications 
p.spc.kappa_0 = 0.8
p.spc.drc_functs_G = [DoseResponse.NEC2neg] 

# save config file
save_config("config/my_config_file.yml", p)
``` 

"""
function save_config(path_to_config::String, p::AbstractParamCollection)::Nothing

    pout = Dict{Any,Any}()

    for sub in fieldnames(typeof(p))
        dictfromsub = convert(Dict, ntfromstruct(eval(:($p.$sub)))) |> Dict{Any,Any}
        deserialize_dict!(dictfromsub)
        pout[sub] = dictfromsub
    end

    YAML.write_file(path_to_config, pout)

    return nothing
end