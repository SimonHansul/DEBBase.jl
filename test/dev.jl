using Parameters
using DifferentialEquations
using Plots

using StaticArrays
using BenchmarkTools


using DEBBase

@benchmark p = simulator(BaseParamCollection())



V1 = [1., 2., 3.]
V2 = [3., 3., 8.]
f(x1, x2) = x1 / x2

@benchmark @. V1 / V2 # 420 ns

@benchmark [x/y for (x,y) in zip(V1,V2)] # 150 ns
