using OrdinaryDiffEq
using Plots


@with_kw mutable struct Params 
    r = .2
end

function du!(du, u, p, t)
    du.x = u.x * p.r
end

condition(u, t, integrator) = 10 < u.x < 1000
effect!(integrator) = integrator.u.a = 1.

u0 = ComponentVector(x = 0.1, a = 0)
tspan = (0,100)
prob = ODEProblem(du!, u0, tspan, Params())
@time sol = solve(prob, callback = DiscreteCallback(condition, effect!));
plot(
    plot(sol.t, [u.x for u in sol.u], yscale = :log10),
    plot(sol.t, [u.a for u in sol.u])
)


du!(
    ::ComponentVector{Float64, Vector{Float64}, Tuple{Axis{(x = 1, a = 2)}}} 
    ::Params
    ::Float64)