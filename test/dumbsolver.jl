#=
A dumb solver for testing purposes, using forward Euler to solve the model. 
This also happens to be necessary for the IBM implementation.-
=#

using Pkg; Pkg.activate("test")
using Revise 
@time using DEBBase
using ComponentArrays

dt = 1/24

#=
Function to update ComponentArray with dumbsolver
=#

function update!(u::Float64, du::Float64, dt::Float64)
    u += du * dt 
end

function update!(u::AbstractVector, du::AbstractVector, dt::Float64)
    update!.(du, u, dt)
end


function update!(u::ComponentArray, du::ComponentArray, dt::Float64)
    for k in keys(u)
        u_k = getproperty(u, k)
        du_k = getproperty(du, k)
        update!(u_k, du_k, dt)
    end
end

Smax_ref = DEBBase.calc_S_max(DEBBaseParams())
p = BaseParamCollection()
deb = p.deb 
deb.Idot_max_rel = 0.1 # 
deb.Idot_max_rel_emb = deb.Idot_max_rel
deb.eta_AS = 0.901

Smax_deb = DEBBase.calc_S_max(deb)
Z = Smax_deb / Smax_ref
deb.X_emb_int = DEBBaseParams().X_emb_int * Z

glb = p.glb
glb.Xdot_in *= Z

u = DEBBase.initialize_statevars(p)
du = copy(u)

u.S
t = 0.
DEBBase.DEB!(du, u, p, t)
u.S
du.S
du.S * dt
update!(u, du, dt)
u.S

u

u.S
du.S

dt = 1/240

