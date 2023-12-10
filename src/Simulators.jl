
function run_model(
    glb::GlobalParams,
    anm::DEBParams
    )

    u0 = [glb.Xdot_in, anm.Xemb_int, anm.Xemb_int * 0.01, 0., 0.]
    tspan = (0, glb.t_max)
    prob = ODEProblem(DEB!, u0, tspan, Params[glb, anm])
    sol = solve(prob)

    return sol
end