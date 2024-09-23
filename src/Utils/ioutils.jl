function extract_colnames(c::R, k::Symbol) where R <: Real
    return k
end

function extract_colnames(c::AbstractVector, k::Symbol)
    return [Symbol("$(k)_$(i)") for i in 1:length(c)]
end

"""
    extract_colnames(u::ComponentVector)
Extract the column names of an output dataframe from an ODE solution object.
"""
function extract_colnames(u::ComponentVector)
    colnames = []
    for k in keys(u)
        push!(colnames, extract_colnames(u[k], k))
    end
    return vcat(colnames...)
end

"""
    extract_colnames(sol::O)::Vector{Symbol} where {O <:ODESolution}

Extract the column names of an output dataframe from an ODE solution object.
"""
extract_colnames(sol::ODESolution)::Vector{Symbol} = vcat([:t], extract_colnames(sol.u[1]))


ymlparse(x::String) = Meta.parse(x) |> eval
ymlparse(x::Vector{String}) = @. Meta.parse(x) |> eval

"""
    load_config(ParamsType::DataType, path_to_config::String)::ParamsType

Load a configuration file and retu
"""
function load_config(ParamsType::DataType, path_to_config::String)::ParamsType
    config = YAML.load_file(path_to_config)
    p = ParamsType()

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