using Parameters
using DifferentialEquations
using Plots

@with_kw mutable struct Params
    mu::Float64 = 0.2
    K::Float64 = 10.0
end

function dN!(du, u, p, t)
    du[1] =  p.mu * u[1] * (1 - (u[1]/p.K))
end

u0 = [1.]
tspan = (0, 100)
prob = ODEProblem(dN!, u0, tspan, Params())
sol = solve(prob)

plot(sol.t, vcat(sol.u...))