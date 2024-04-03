
using BenchmarkTools
using DEBBase

out = DEBBase.simulator(DEBParamCollection())
@benchmark out = DEBBase.simulator(
        DEBParamCollection(
            glb = GlobalParams(Xdot_in = 4800., t_max = 21.), 
            deb = SpeciesParams(K_X = 12e3))
        )