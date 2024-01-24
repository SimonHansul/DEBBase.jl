
"""
Sigmoid switch function. 
`y_left` and `y_right` are the function values left and right of the threshold `x_thr`.
$(TYPEDSIGNATURES)
"""
@inline function sig(
    x::Float64, 
    x_thr::Float64,
    y_left::Float64, 
    y_right::Float64; 
    β::Float64 = 1e6
    )
    return 1 / (1 + exp(-β*(x - x_thr))) * (y_right - y_left) + y_left
end

@inline function functional_response(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )
    let X_V = u.X_p / p.glb.V_patch
        return X_V / (X_V + p.deb.K_X)
    end
end

"""
Calculate ingestion rate. 
Embryos (life_stage ≈ 1.) take up resources from the vitellus X_emb. 
Juveniles and adults (life_stage > 1) feed on the external resource X_p.
$(TYPEDSIGNATURES)
"""
@inline function Idot!(
    du::ComponentArray,
    u::ComponentArray, 
    p::AbstractParamCollection, 
    t::Real
    )

    du.I_emb = sig(
        u.X_emb, # uptake from vitellus depends on mass of vitellus
        0., # the switch occurs when vitellus is used up 
        0., # when the vitellus is used up, there is no uptake
        u.S^(2/3) * p.deb.Idot_max_rel # when the vitellus is not used up, uptake from vitellus occurs
        )
    du.I_p = sig(
        u.X_emb, # ingestion from external resource depends on mass of vitellus
        0., # the switch occurs when the vitellus is used up  
        functional_response(du, u, p, t) * p.deb.Idot_max_rel * u.S^(2/3), # when the vitellus is used up, ingestion from the external resource occurs
        0. # while there is still vitellus left, there is no uptake from the external resource
        )
    du.I = du.I_emb + du.I_p
end

"""
Assimilation rate
$(TYPEDSIGNATURES)
"""
@inline function Adot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )
    du.A = du.I * p.deb.eta_IA * u.y_A
end

"""
Somatic maintenance rate
$(TYPEDSIGNATURES)
"""
@inline function Mdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real)
    du.M = u.S * p.deb.k_M * u.y_M
end

"""
Maturity maintenance rate
$(TYPEDSIGNATURES)
"""
@inline function Jdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )
    du.J = u.H * p.deb.k_J
end

"""
Positive somatic growth
$(TYPEDSIGNATURES)
"""
@inline function Sdot_positive(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )
    return p.deb.eta_AS * u.y_G * (p.deb.kappa * du.A - du.M)
end

"""
Negative somatic growth
($(TYPEDSIGNATURES))
"""
@inline function Sdot_negative(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )
    return -(du.M / p.deb.eta_SA - p.deb.kappa * du.A)
end

"""
Somatic growth rate, including application of shrinking equation.
$(TYPEDSIGNATURES)
"""
@inline function Sdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )
    du.S = (p.deb.kappa * u.A >= du.M) * Sdot_positive(du, u, p, t) + (p.deb.kappa * u.A < du.M) * Sdot_negative(du, u, p, t)
end

"""
Maturation rate. 
Maturity is dissipated energy and can therefore not be burned to cover maintenance costs. <br>
Currently there are no consequences for an organism not covering maturity maintenance. 
If this turns out to be an issue, we might consider to add a damage pool D_H,
where the amount of maturity maintenance that could not be covered is accumulated. 
This might then lead to a fitness penalty depending on D_H, for example in the form of additional 
mortality or embryonic hazard (TBD). 
$(TYPEDSIGNATURES)
"""
@inline function Hdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )
    du.H = sig(
        u.H, # maturation depends on maturity
        p.deb.H_p, # switch occurs at maturity at puberty H_p
        # TODO: replace max(0, .) with sigmoid function
        max(0., ((1 - p.deb.kappa) * du.A) - du.J), # maturation for embryos and juveniles
        0., # maturation for adults
    )
end

"""
Reproduction rate.
$(TYPEDSIGNATURES)
"""
@inline function Rdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )
    du.R = sig(
        u.H, # reproduction depends on maturity
        p.deb.H_p, # switch occurs at maturity at puberty H_p
        0., # reproduction for embryos and juveniles
        u.y_R * p.deb.eta_AR * (1 - p.deb.kappa) * du.A - du.J # reproduction for adults
    )
end

"""
Calculation of the total dissipation flux, equal to the sum of maintenance costs and overheads paid for assimilation, mobilization, growth and reproduction.
$(TYPEDSIGNATURES)
"""
function Qdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real; 
    Idot::Float64, 
    Mdot::Float64, 
    Jdot::Float64,
    )
    # dissipation fluxes for the individual processes
    let Qdot_A, Qdot_S, Qdot_C, Qdot_R
        Qdot_A = Idot * (1 - deb.eta_IA)
        Qdot_S = du.S >= 0 ? du.S * (1 - p.deb.eta_AS)/p.deb.eta_AS : du.S * (p.deb.eta_SA - 1)
        Qdot_R = du.R * (1 - p.deb.eta_AR)/p.deb.eta_AR
        du.Q =  Qdot_A + Qdot_S + Qdot_C + Qdot_R + Mdot + Jdot + Hdot
    end
end


"""
TK for DEB-TKTD model, including effect of surface area to volume ratio and dilution by growth. 
If `D` and is given as a Vector, TK is only stressor-specific but not PMoA-specific. 
$(TYPEDSIGNATURES)
"""
@inline function Ddot!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    ) 


    # TODO: move calculation of L_S_max out so it is only calculated once, not at every step
    S_max = calc_S_max(p.deb) # maximum structural length (g^(1/3))
    AV_min = S_max^(2/3) / S_max # minimum area-to-volume (AV) ratio is reached at maximum size
    AV = u.S^(2/3) / u.S # current AV ratio
    AV_rel = AV / AV_min # correction factor is AV ratio relative to minimum AV ratio 
    for z in eachindex(u.C_W)
        du.D_G[z] = sig(u.X_emb, 0., p.deb.k_D_G[z] * AV_rel * (u.C_W[z] - u.D_G[z]) - u.D_G[z] * (du.S / u.S), 0.)
        du.D_M[z] = sig(u.X_emb, 0., p.deb.k_D_M[z] * AV_rel * (u.C_W[z] - u.D_M[z]) - u.D_M[z] * (du.S / u.S), 0.)
        du.D_A[z] = sig(u.X_emb, 0., p.deb.k_D_A[z] * AV_rel * (u.C_W[z] - u.D_A[z]) - u.D_A[z] * (du.S / u.S), 0.)
        du.D_R[z] = sig(u.X_emb, 0., p.deb.k_D_R[z] * AV_rel * (u.C_W[z] - u.D_R[z]) - u.D_R[z] * (du.S / u.S), 0.)
        du.D_h[z] = sig(u.X_emb, 0., p.deb.k_D_h[z] * AV_rel * (u.C_W[z] - u.D_h[z]) - u.D_h[z] * (du.S / u.S), 0.)
    end
end

@inline function C_Wdot!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )
    du.C_W = zeros(length(u.C_W)) # constant exposure : derivative is 0
end

@inline function X_pdot!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )
    du.X_p = p.glb.Xdot_in - du.I_p
end

@inline function X_embdot!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )
    du.X_emb = -du.I_emb
end

@inline function y!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )
    u.y_G = prod([p.deb.drc_functs_G[z](u.D_G[z], p.deb.drc_params_G[z]) for z in 1:length(u.C_W)])
    u.y_M = prod([p.deb.drc_functs_M[z](u.D_M[z], p.deb.drc_params_M[z]) for z in 1:length(u.C_W)])
    u.y_A = prod([p.deb.drc_functs_A[z](u.D_A[z], p.deb.drc_params_A[z]) for z in 1:length(u.C_W)])
    u.y_R = prod([p.deb.drc_functs_R[z](u.D_R[z], p.deb.drc_params_R[z]) for z in 1:length(u.C_W)])
    u.h_z = sum([p.deb.drc_functs_h[z](u.D_h[z], p.deb.drc_params_h[z]) for z in 1:length(u.C_W)])
end


"""
Definition of reserveless DEB derivatives. 
$(TYPEDSIGNATURES)
"""
function DEB!(du, u, p, t)
    #### boilerplate
    u.S = sig(u.S, 0., 0., u.S)
    
    #### stressor responses
    y!(du, u, p, t)

    #### auxiliary state variables (record cumulative values)
    Idot!(du, u, p, t)
    Adot!(du, u, p, t) 
    Mdot!(du, u, p, t) 
    Jdot!(du, u, p, t)

    #### major state variables
    Sdot!(du, u, p, t) # structure
    Hdot!(du, u, p, t) # maturity 
    Rdot!(du, u, p, t) # reproduction buffer
    X_pdot!(du, u, p, t) # resource abundance
    X_embdot!(du, u, p, t) # vitellus
    Ddot!(du, u, p, t) # damage
    C_Wdot!(du, u, p, t) # external stressor concentration  
end
