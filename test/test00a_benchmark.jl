
using Pkg; Pkg.activate("test")
using BenchmarkTools
using DEBBase
using OrdinaryDiffEq
using Test

@benchmark out = DEBBase.simulator(DEBParamCollection(), alg = Tsit5()) 
@test true # this just has to run without an error