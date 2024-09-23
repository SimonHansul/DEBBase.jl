using Pkg; Pkg.activate("test")

using YAML
using Revise
using DEBBase.DEBODE
using DEBBase.DoseResponse
using DEBBase.Utils
using Distributions
using Test

p = Params()
p.spc.Z = LogNormal(1, 1)

using NamedTupleTools
pout = Dict{Any,Any}()

deserialize(x::Any) = x # most types don't need extra treatment

function deserialize(d::Distribution)
    disttype = "$(typeof(d))"
    distparams = fieldnames(typeof(d))
    deserialized_dist_params = [deserialize(getproperty(d, p)) for p in distparams]
    return "$disttype($(["$(x), " for x in deserialized_dist_params]...))"
end

d = Normal(1, 1) #Truncated(LogNormal(1, 0.1), 0, Inf)
@test deserialize(d) == "Normal{Float64}(1.0, 1.0, )"

function deserialize(d::Truncated)
    untruncated = deserialize(d.untruncated)
    return "Truncated($untruncated, $(d.lower), $(d.upper))"
end

d = Truncated(Normal(1, 1), 0, Inf)
deserialize(d)

d.upper

function deserialize_dict!(dict)::Nothing

    for key in keys(dict)
        dict[key] = deserialize_val(dict[key])
    end

    return nothing 
end


for sub in fieldnames(typeof(p))
    dictfromsub = convert(Dict, ntfromstruct(eval(:(p.$sub))))
    deserialize_dict!(dictfromsub)
    pout[sub] = dictfromsub
end

pout

pout

YAML.write_file("test/config/save_config_example.yml", pout)

# saving the file works, but s
load_config(Params, "test/config/save_config_example.yml")


NamedTuple(x = 1, 2 = 3)
using Parameters

@with_kw mutable struct mstr
    x::NamedTuple
end

mstr((a = 2, b = 3))