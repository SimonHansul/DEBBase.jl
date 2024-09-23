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


p.spc.drc_functs_G