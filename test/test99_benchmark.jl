using DEBBase.DEBODE

b = @benchmark yhat = simulator(DEBParamCollection()) 
@info("Median benchmark at $(median(b.times)/1e6) ms for default parameters")
@test median(b.times) < 5e6 # computation time should be 5ms for the default parameters

yhat = simulator(DEBParamCollection())
yhat

@df yhat plot(
    plot(:t, :X_emb), 
    plot(:t, :S),
    plot(:t, :H),
    plot(:t, [:embryo, :juvenile, :adult], leg = true, label = ["embryo" "juvenile" "adult"])
)