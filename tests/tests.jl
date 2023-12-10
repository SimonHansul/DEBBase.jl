using Revise
@time using DEBBase
using ProfileView

default( # plotting defaults
    legend = false,
    linewidth = 1.5,
    labelfontsize = 10,
    legendtitlefontsize = 10,
    fillalpha = .2,
    foreground_color_legend = nothing
    )

glb = DEBBase.GlobalParams()
glb.Xdot_in = 1e10

anm = DEBBase.DEBParams(k_J = 0.)



#FIXME : model implementation dies not respond properly to changes in parameters ??
#TODO  : check whether Pararam types cause considerable overhead
@time sol = DEBBase.run_model(glb, anm);
plot((hcat(out.u...)'[:,3]))


DEBBase.Params[glb, anm]