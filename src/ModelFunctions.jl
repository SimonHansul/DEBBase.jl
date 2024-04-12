"""
Clip negative values.
"""
function clipneg(x::Float64)
    return sig(x, 0., 0., x)
end

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
    beta::Float64 = 1e16
    )
    return 1 / (1 + exp(-beta*(x - x_thr))) * (y_right - y_left) + y_left
end

@inline function functional_response(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    ) 
    let X_V = u.X_p / p.glb.V_patch # convert food abundance to concentration
        return X_V / (X_V + p.agn.K_X) # calculate type II functional response
    end
end

"""
Calculate ingestion rate. 
Embryos (X_emb <= 0) take up resources from the vitellus X_emb. 
Juveniles and adults (X_emb > 0) feed on the external resource X_pcmn.
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
        (Complex(u.S)^(2/3)).re * p.agn.Idot_max_rel; # when the vitellus is not used up, uptake from vitellus occurs
        beta = 1e20 # for switches around 0, we need very high beta values
        )

    du.I_p = sig(
        u.X_emb, # ingestion from external resource depends on mass of vitellus
        0., # the switch occurs when the vitellus is used up  
        functional_response(du, u, p, t) * p.agn.Idot_max_rel * (Complex(u.S)^(2/3)).re, # when the vitellus is used up, ingestion from the external resource occurs
        0.; # while there is still vitellus left, there is no uptake from the external resource
        beta = 1e20 # again we have a switch around 0, requiring very high beta
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
    t::Real) 
    du.A = du.I * p.spc.eta_IA * u.y_A
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
    du.M = u.S * p.spc.k_M * u.y_M
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
    du.J = u.H * p.spc.k_J * u.y_M
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
    return p.spc.eta_AS * u.y_G * (p.spc.kappa * du.A - du.M)
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
    return -(du.M / p.spc.eta_SA - p.spc.kappa * du.A)
end

function Sdot(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    ) 
    return sig(
        p.spc.kappa * du.A, # growth depends on maintenance coverage
        du.M, # switch occurs based on maintenance costs
        Sdot_negative(du, u, p, t), # left of the threshold == maintenance costs cannot be covered == negative growth
        Sdot_positive(du, u, p, t) # right of the threshold == maintenance costs can be covered == positive growth
    )
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
    du.S = Sdot(du, u, p, t)
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
        p.agn.H_p, # switch occurs at maturity at puberty H_p
        clipneg(((1 - p.spc.kappa) * du.A) - du.J), # maturation for embryos and juveniles
        0., # maturation for adults
    )
end

"""
Update the current estimate of H_b. 
The current estimate of maturity at birth is equal to current maturity for embryos, 
and will be fixed to the current value upon completion of embryonic development. \n
This way, we obtain H_b as an implied trait and can use it later (for example in `abj()`).
"""
@inline function H_bdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    ) 
    du.H_b = DEBBase.sig(
        u.X_emb, # estimate depends on embryonic state
        0., # switch occurs when vitellus is gone
        0., # post-embryonic: H_b stays fixed
        du.H # embryonic: H_b tracks H
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
        p.agn.H_p, # switch occurs at maturity at puberty H_p
        0., # reproduction for embryos and juveniles
        clipneg(u.y_R * p.spc.eta_AR * ((1 - p.spc.kappa) * du.A - du.J)) # reproduction for adults
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
    t::Real
    ) 
    # dissipation fluxes for the individual processes
    let Qdot_A, Qdot_S, Qdot_C, Qdot_R
        Qdot_A = Idot * (1 - p.spc.eta_IA)
        Qdot_S = du.S >= 0 ? du.S * (1 - p.spc.eta_AS)/p.spc.eta_AS : du.S * (p.spc.eta_SA - 1)
        Qdot_R = du.R * (1 - p.spc.eta_AR)/p.spc.eta_AR
        du.Q =  Qdot_A + Qdot_S + Qdot_C + Qdot_R + Mdot + Jdot + Hdot
    end
end


"""
TK for spc-TKTD model, including effect of surface area to volume ratio and dilution by growth. 
If `D` and is given as a Vector, TK is only stressor-specific but not PMoA-specific. 
$(TYPEDSIGNATURES)
"""
@inline function Ddot!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    ) 

    for z in eachindex(u.C_W)
        # the sigmoid function causes Ddot to be 0 for embryos (assumption of internal eggs which are not exposed to external stressor )
        @inbounds du.D_G[z] = sig(u.X_emb, 0., p.spc.k_D_G[z] * (u.C_W[z] - u.D_G[z]), 0.)
        @inbounds du.D_M[z] = sig(u.X_emb, 0., p.spc.k_D_M[z] * (u.C_W[z] - u.D_M[z]), 0.)
        @inbounds du.D_A[z] = sig(u.X_emb, 0., p.spc.k_D_A[z] * (u.C_W[z] - u.D_A[z]), 0.)
        @inbounds du.D_R[z] = sig(u.X_emb, 0., p.spc.k_D_R[z] * (u.C_W[z] - u.D_R[z]), 0.)
        @inbounds du.D_h[z] = sig(u.X_emb, 0., p.spc.k_D_h[z] * (u.C_W[z] - u.D_h[z]), 0.)
        
        #@inbounds du.D_G[z] = sig(u.X_emb, 0., p.spc.k_D_G[z] * (calc_SL_max(p.spc) / (Complex(u.S)^(1/3)).re) * (u.C_W[z] - u.D_G[z]) - u.D_G[z] * (du.S / u.S), 0.)
        #@inbounds du.D_M[z] = sig(u.X_emb, 0., p.spc.k_D_M[z] * (calc_SL_max(p.spc) / (Complex(u.S)^(1/3)).re) * (u.C_W[z] - u.D_M[z]) - u.D_M[z] * (du.S / u.S), 0.)
        #@inbounds du.D_A[z] = sig(u.X_emb, 0., p.spc.k_D_A[z] * (calc_SL_max(p.spc) / (Complex(u.S)^(1/3)).re) * (u.C_W[z] - u.D_A[z]) - u.D_A[z] * (du.S / u.S), 0.)
        #@inbounds du.D_R[z] = sig(u.X_emb, 0., p.spc.k_D_R[z] * (calc_SL_max(p.spc) / (Complex(u.S)^(1/3)).re) * (u.C_W[z] - u.D_R[z]) - u.D_R[z] * (du.S / u.S), 0.)
        #
        #@inbounds du.D_h[z] = sig(u.X_emb, 0., p.spc.k_D_h[z] * (calc_SL_max(p.spc) / (Complex(u.S)^(1/3)).re) * (u.C_W[z] - u.D_h[z]) - u.D_h[z] * (du.S / u.S), 0.)
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

    @inbounds u.y_G = prod([p.spc.drc_functs_G[z](u.D_G[z], (p.spc.e_G[z], p.spc.b_G[z])) for z in 1:length(u.C_W)]) # combined relative responses for sublethal effects
    @inbounds u.y_M = prod([p.spc.drc_functs_M[z](u.D_M[z], (p.spc.e_M[z], p.spc.b_M[z])) for z in 1:length(u.C_W)])
    @inbounds u.y_A = prod([p.spc.drc_functs_A[z](u.D_A[z], (p.spc.e_A[z], p.spc.b_A[z])) for z in 1:length(u.C_W)])
    @inbounds u.y_R = prod([p.spc.drc_functs_R[z](u.D_R[z], (p.spc.e_R[z], p.spc.b_R[z])) for z in 1:length(u.C_W)])

    @inbounds u.h_z = sum([p.spc.drc_functs_h[z](u.D_h[z], (p.spc.e_h[z], p.spc.b_h[z])) for z in 1:length(u.C_W)]) # hazard rate
end

"""
Metabolic acceleration from birth to maturity. 
We assume that some baseline parameter `p` has value `p_b` at birth and `p_j` at metamorphosis.
Between birth and metamorphosis, the current value of `p` is the maturity-weighted mean of `p_b` and `p_j`.
"""
function Hbj(H::Float64, X_emb::Float64, H_b::Float64, H_j::Float64, p_b::Float64, p_j::Float64)::Float64
    w_b = (H_j - H) / (H_j - H_b) # weight for p_b
    w_j = 1 - w_b # weight for p_j
    p_bj = mean([p_b, p_j], Weights([w_b, w_j])) # p_bj, i.e. value between birth and maturity
    
    p = DEBBase.sig( # post-metamorphosis: value stays constant at p_j
        H,
        H_j,
        DEBBase.sig( # embryonic: value stays constant at p_b
            X_emb, 
            0., 
            p_bj, 
            p_b),
        p_j
    )
    return p
end

"""
Definition of base model system. 
$(TYPEDSIGNATURES)
"""
function DEB!(du, u, p, t)

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
    H_bdot!(du, u, p, t) # estimate of maturity at birth
    Rdot!(du, u, p, t) # reproduction buffer
    X_pdot!(du, u, p, t) # resource abundance
    X_embdot!(du, u, p, t) # vitellus
    Ddot!(du, u, p, t) # damage
    C_Wdot!(du, u, p, t) # external stressor concentration  
end