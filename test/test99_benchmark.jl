using BenchmarkTools
using ProfileView

using DEBBase.DEBODE

b = @benchmark yhat = simulator(DEBParamCollection()) 
@info("Median benchmark at $(median(b.times)/1e6) ms for default parameters")
@test median(b.times) < 10e6 # computation time should be below 5ms for the default parameters

