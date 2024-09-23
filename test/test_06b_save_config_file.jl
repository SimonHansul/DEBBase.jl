using Pkg; Pkg.activate("test")

using YAML
using Revise
using DEBBase.DEBODE
using DEBBase.DoseResponse
using DEBBase.Utils
using Distributions
using Test


@testset begin
    @info "Trying to load the example configuration file "
    p = load_config(Params, "test/config/config_example.yml")
    @info "Checking return type"
    @test p isa Params
    @info "Checking values of the returned object"
    @test p.spc.Z == Truncated(Normal(1.33, 0.133), 0, Inf)
    @test p.spc.drc_functs_G == [DoseResponse.NEC2neg]
end


p.spc.Z = LogNormal(1, 1)

using NamedTupleTools

sub = :glb




pout = Dict{Any,Any}()

for sub in fieldnames(typeof(p))
    pout[sub] = convert(Dict, ntfromstruct(eval(:(p.$sub))))
end

pout

YAML.write_file("test/config/save_config_example.yml", pout)



load_config(Params, "test/config/save_config_example.yml")


NamedTuple(x = 1, 2 = 3)
using Parameters

@with_kw mutable struct mstr
    x::NamedTuple
end

mstr((a = 2, b = 3))