using Infiltrator
using Revise
@time using DEBBase

#FIXME : model implementation dies not respond properly to changes in parameters ??
#TODO  : check whether Pararam types cause considerable overhead

#### 

glb = GlobalBaseParams()
anm = DEBBaseParams()
sol = DEBBase.run_model(glb, anm)

using Plots
plot(sol.t, hcat(sol.u...)'[:,3])


DEBBase.calc_S_max(anm)