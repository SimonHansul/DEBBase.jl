
"""
Run the DEBBase Model. 
$(TYPEDSIGNATURES)
"""
function run_model(
    glb::GlobalBaseParams,
    anm::DEBBaseParams
    )

    u0 = [glb.Xdot_in, anm.X_emb_int, anm.X_emb_int * 0.01, 0., 0.]
    tspan = (0, glb.t_max)
    prob = ODEProblem(DEB!, u0, tspan, AbstractParams[glb, anm])
    sol = solve(prob, reltol = 1e-6, abstol = 1e-10)

    return sol
end