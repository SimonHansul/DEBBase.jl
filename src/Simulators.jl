
"""
Run the DEBBase Model. 
$(TYPEDSIGNATURES)
"""
function run_model(
    glb::AbstractParams,
    deb::AbstractParams
    )

    u0 = [glb.Xdot_in, deb.X_emb_int, deb.X_emb_int * 0.01, 0., 0.]
    tspan = (0, glb.t_max)
    prob = ODEProblem(DEB!, u0, tspan, (glb = glb, deb = deb))
    sol = solve(prob, reltol = 1e-6, abstol = 1e-10)
    simout = DataFrame(hcat(sol.t, hcat(sol.u...)'), [:t, :X, :X_emb, :S, :H, :R])

    return simout
end