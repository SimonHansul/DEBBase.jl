using Pkg; Pkg.activate("tests")
using Revise
using DEBBase
using BenchmarkTools

BaseParamCollection().deb

# 3.4 MiB memeory, 9346 allocs, 7.5 ms
# 3 MiB, 81050 allocs, 6.6 ms
@benchmark simulator(BaseParamCollection())

using DEBBase
out = simulator(BaseParamCollection())
