
"""
Run any DEBBase-compatible model. 
$(TYPEDSIGNATURES)
"""
function simulator(
    glb::AbstractParams,
    deb::AbstractParams
    )
    
    u = ComponentArray(
        X_p = glb.Xdot_in, 
        X_emb = deb.X_emb_int, 
        S = deb.X_emb_int * 0.01, 
        H = 0., 
        R = 0.,
        D = zeros(length(glb.C_W))
    )
    tspan = (0, glb.t_max)
    prob = ODEProblem(DEB!, u, tspan, (glb = glb, deb = deb))
    sol = solve(prob, reltol = 1e-6, abstol = 1e-10)
    simout = DataFrame(hcat(sol.t, hcat(sol.u...)'), [:t, :X, :X_emb, :S, :H, :R])

    return simout
end