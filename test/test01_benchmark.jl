using Pkg; Pkg.activate("test")
using BenchmarkTools
using OrdinaryDiffEq
using Test
using Revise
using DEBBase
using DEBBase.DEBODE

b = @benchmark out = simulator(DEBParamCollection()) 
@test median(b.times) < 10e6 # computation time should be (clearly) below 10ms for the default parameters