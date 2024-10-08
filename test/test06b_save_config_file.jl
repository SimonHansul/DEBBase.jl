using Pkg; Pkg.activate("test")

using Revise
using DEBBase.DEBODE
using DEBBase.DoseResponse
using DEBBase.Utils
using Distributions
using Test

using DataStructures
using NamedTupleTools

@info "Deserializing a distribution"
d = LogNormal(1, 0.1)
@test Utils.deserialize(d) == "LogNormal{Float64}(1.0, 0.1, )"

@info "Deserializing a truncated distribution"
d = Truncated(Normal(1, 1), 0, Inf)
@test Utils.deserialize(d) == "Truncated(Normal{Float64}(1.0, 1.0, ), 0.0, Inf)"

@info "Loading a config file and checking values"
pload = params_from_config(Params, "test/config/param_config_save_example.yml")
@test pload.spc.Z == LogNormal(1, 1) 


# what if we define a new DRC and try to save/load it?
# does currently not work, TBD
#newdrc(x, p) = x*p[1] + p[2]
#
#p = Params()
#p.spc.drc_functs_M = [newdrc]
#save_config("./test/config/config_with_custom_drc.yml", p)
#
#params_from_config(Params, "./test/config/config_with_custom_drc.yml")