using BenchmarkTools
using DEBBase.DEBODE

# trying 

yhat = DEBODE.simulator(Params())

b = @benchmark yhat = DEBODE.simulator(Params()) 
@info("Median benchmark at $(median(b.times)/1e6) ms for default parameters")
@test median(b.times) < 10e6 # computation time should be below 5ms for the default parameters

DEBODE.treplicates

@benchmark DEBODE.treplicates(DEBODE.simulator, Params(), 20)

#using ProfileView
#VSCodeServer.@profview [DEBODE.simulator(Params())  for _ in 1:100]

using DEBBase.AgentBased
@benchmark yhat = AgentBased.simulator(Params())