
"""
Run any DEBBase-compatible model. 
$(TYPEDSIGNATURES)
"""
function simulator(
    glb::AbstractParams,
    deb::AbstractParams
    )
    
    @assert size(deb.k_D)[1] == length(glb.C_W) "Size of k_D ($(size(deb.k_D))) does not match number of stressors indicated in glb.C_W ($(length(glb.C_W)))."
    @assert size(deb.drc_functs)[1] >= length(glb.C_W) "Size of drc_functs ($(size(deb.drc_functs))) does not match number of stressors indicated in glb.C_W ($(length(glb.C_W)))."
    @assert size(deb.drc_params)[1] >= length(glb.C_W) "Size of drc_params ($(size(deb.drc_params))) does not match number of stressors indicated in glb.C_W ($(length(glb.C_W)))."

    u = ComponentArray( # initial states
        X_p = glb.Xdot_in, # initial resource abundance equal to influx rate
        X_emb = deb.X_emb_int, # initial mass of vitellus
        S = deb.X_emb_int * 0.01, # initial structure is a small fraction of initial reserve // mass of vitellus
        H = 0., # maturity
        R = 0., # reproduction buffer
        D = zeros(size(deb.k_D)),
        C_W = glb.C_W
    )

    tspan = (0, glb.t_max) # define the timespan
    prob = ODEProblem(DEB!, u, tspan, (glb = glb, deb = deb)) # define the problem
    
    # TODO: flatten D-Matrix into multiple columns in output dataframe
    sol = solve(prob, reltol = 1e-6, abstol = 1e-10) # get solution to the IVP

    simout = DataFrame(
        hcat(sol.t, hcat(sol.u...)'), # ODE output converted to wide matrix
        vcat(  # column names
        #:# TODO: assign column names for different PMoAs dynamically
        [:t, :X, :X_emb, :S, :H, :R],
        [:D_1_1, :D_2_1],
        [:C_W_1]
        ))

    return simout
end