using BenchmarkTools


yhat = simulator(Params())

b = @benchmark yhat = simulator(Params()) 
@info("Median benchmark at $(median(b.times)/1e6) ms for default parameters")
@test median(b.times) < 10e6 # computation time should be below 5ms for the default parameters

using ProfileView
VSCodeServer.@profview [simulator(Params())  for _ in 1:100]