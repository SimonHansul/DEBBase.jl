using BenchmarkTools
using ProfileView

using DEBBase.DEBODE

VSCodeServer.@profview [simulator(DEBParamCollection()) for _ in 1:100]

b = @benchmark yhat = simulator(DEBParamCollection(), alg = Tsit5()) 
@info("Median benchmark at $(median(b.times)/1e6) ms for default parameters")
@test median(b.times) < 10e6 # computation time should be below 5ms for the default parameters

yhat = simulator(DEBParamCollection())
yhat

@df yhat plot(
    plot(:t, :X_emb), 
    plot(:t, :S),
    plot(:t, :H),
    plot(:t, [:embryo, :juvenile, :adult], leg = true, label = ["embryo" "juvenile" "adult"])
)

