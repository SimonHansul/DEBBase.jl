using BenchmarkTools
using DEBBase.DEBODE
using OrdinaryDiffEq
ROCK2

@time yhat = DEBODE.simulator(Params(), alg = Tsit5());
@time yhat = DEBODE.simulator(Params(), alg = Tsit5(), reltol = 1e-4);

@benchmark DEBODE.simulator(Params(), alg = Tsit5())

@benchmark DEBODE.simulator(p, reltol = 1e-2)


begin
    params.agn = DEBODE.AgentParamType(params.spc) # initialize agent parameters incl. individual variability
    callbacks = DEBODE.lifestage_callbacks()
    
    u = DEBODE.initialize_statevars(params)
    prob = ODEProblem(model, u, (0, params.glb.t_max), params) # define the problem
    sol = solve(prob; callback = callbacks, saveat = saveat, reltol = reltol, kwargs...) # get solution to the IVP
    simout = DEBODE.sol_to_df(sol) # convert solution to dataframe
    
    b = @benchmark yhat = DEBODE.simulator(Params()) 
    @info("Median benchmark at $(median(b.times)/1e6) ms for default parameters")
    @test median(b.times) < 10e6 # computation time should be below 5ms for the default parameters    
end

DEBODE.treplicates

@benchmark DEBODE.treplicates(DEBODE.simulator, Params(), 20)

#using ProfileView
#VSCodeServer.@profview [DEBODE.simulator(Params())  for _ in 1:100]

using DEBBase.DEBABM
@benchmark yhat = DEBABM.simulator(Params())