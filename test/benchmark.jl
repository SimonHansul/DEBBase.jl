
using BenchmarkTools
using DEBBase

out = Simulator(BaseParamCollection())
@benchmark out = simulator(
        BaseParamCollection(
            glb = GlobalBaseParams(Xdot_in = 4800., t_max = 21.), 
            deb = DEBBaseParams(K_X = 12e3))
        )