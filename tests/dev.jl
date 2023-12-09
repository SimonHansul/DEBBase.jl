using DifferntialEquations

function du!(du, u, p, t)
    return u[1] * p["mu"]
end

u0 = (1)
tspan = (0,3)
p = Dict(:mu => 0.2)

prob = ODEProblem(du!, u0, tspan, p) # define the initial value problem
sol = solve(prob, alg, reltol = reltol, saveat = glb[:saveat]; kwargs...) # solve the IVP
