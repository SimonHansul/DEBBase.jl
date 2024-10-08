using Pkg; Pkg.activate("test")

using Revise
using DEBBase.DEBODE
using DEBBase.DoseResponse
using DEBBase.Utils
using Distributions
using Test

@testset begin
    @info "Trying to load the example parameter configuration file "
    p = params_from_config(Params, "test/config/param_config_example.yml")
    @info "Checking return type"
    @test p isa Params
    @info "Checking values of the returned object"
    @test p.spc.Z == Truncated(Normal(1.33, 0.133), 0, Inf)
    @test p.spc.drc_functs_G == [DoseResponse.NEC2neg]
    @test p.spc.e_G == [1.0]
    @test p.spc.b_G == [2.0]
end
