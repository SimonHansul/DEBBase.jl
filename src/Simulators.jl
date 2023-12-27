
"""
Run any DEBBase-compatible model. 
$(TYPEDSIGNATURES)
"""
function simulator(
    glb::AbstractParams,
    deb::AbstractParams
    )
    
    @assert size(deb.k_D)[2] == length(glb.C_W) "Number of rows in deb.k_D ($(size(deb.k_D)[2])) has to match number of stressors indicated in glb.C_W ($(length(glb.C_W)))."

    u = ComponentArray( # initial states
        X_p = glb.Xdot_in, # initial resource abundance equal to influx rate
        X_emb = deb.X_emb_int, # initial mass of vitellus
        S = deb.X_emb_int * 0.01, # initial structure is a small fraction of initial reserve // mass of vitellus
        H = 0., # maturity
        R = 0., # reproduction buffer
        D = zeros(size(deb.k_D))
    )

    tspan = (0, glb.t_max) # define the timespan
    prob = ODEProblem(DEB!, u, tspan, (glb = glb, deb = deb)) # define the problem
    
    # TODO: flatten D-Matrix into multiple columns in output dataframe
    sol = solve(prob, reltol = 1e-6, abstol = 1e-10) # get solution to the IVP

    simout = DataFrame(
        hcat(sol.t, hcat(sol.u...)'), # ODE output converted to wide matrix
        [:t, :X, :X_emb, :S, :H, :R, :D_G, :D_M, :D_A, :D_R, :D_H] # column names
        )
    
    return simout
end