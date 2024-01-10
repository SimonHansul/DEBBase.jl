
"""
Sigmoid switch function. 
`y_left` and `y_right` are the function values left and right of the threshold `x_thr`.
$(TYPEDSIGNATURES)
"""
function sig(
    x::Float64, 
    x_thr::Float64,
    y_left::Float64, 
    y_right::Float64; 
    β::Float64 = 100.
    )
    return 1 / (1 + exp(-β*(x - x_thr))) * (y_right - y_left) + y_left
end


"""
Determine the current life stage. 
$(TYPEDSIGNATURES)
"""
@inline function determine_life_stage!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )
    if u.H >= p.deb.H_p
        u.life_stage = 3 # adult
    elseif u.X_emb <= 0
        u.life_stage = 2 # juvenile
    else
        u.life_stage = 1 # embryo
    end
end


"""
$(TYPEDSIGNATURES)
"""
@inline function embryo(life_stage::Float64)
    return life_stage == 1
end

"""
$(TYPEDSIGNATURES)
"""
@inline function juvenile(life_stage::Float64)
    return life_stage == 2
end

"""
$(TYPEDSIGNATURES)
"""
@inline function adult(life_stage::Float64)
    return life_stage == 3
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
Calculate ingestion rate for a single resource.
$(TYPEDSIGNATURES)
"""
@inline function Idot!(
    du::ComponentArray,
    u::ComponentArray, 
    p::AbstractParamCollection, 
    t::Real
    )
    
    if u.life_stage > 1 # juveniles and adults feed from external resource
        du.I = functional_response(du, u, p, t) * p.deb.Idot_max_rel_emb * u.S^(2/3) 
    else # embryos feed from the vitellus
        du.I = u.S^(2/3) * p.deb.Idot_max_rel
    end
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
    du.M = max(0, u.S * p.deb.k_M * u.y_M)
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
Maturity is dissipated energy and can therefore not be burned to cover maintenance costs. 
For simplicity, we currently assume that maturity maintenance will not be covered of the 1-kappa flux is insufficient.
$(TYPEDSIGNATURES)
"""
@inline function Hdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )
    if !adult(u.life_stage)
        du.H =  max(0, ((1 - p.deb.kappa) * du.A) - du.J)
    else
        du.H = 0.0
    end
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
    if adult(u.life_stage)
        du.R = u.y_R * p.deb.eta_AR * (1 - p.deb.kappa) * du.A - du.J 
    else
        du.R =  0.0
    end
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

    let SL = u.S^(1/3) / u.S, # strucutral length (g^(1/3))
        # TODO: move calculation of L_S_max out so it is only calculated once, not at every step
        SL_max = calc_SL_max(p.deb) # maximum structural length (g^(1/3))
        @. du.D_G = p.deb.k_D_G * (SL_max / SL) * (u.C_W - u.D_G) - u.D_G * (du.S / u.S)
        @. du.D_M = p.deb.k_D_M * (SL_max / SL) * (u.C_W - u.D_M) - u.D_M * (du.S / u.S)
        @. du.D_A = p.deb.k_D_A * (SL_max / SL) * (u.C_W - u.D_A) - u.D_A * (du.S / u.S)
        @. du.D_R = p.deb.k_D_R * (SL_max / SL) * (u.C_W - u.D_R) - u.D_R * (du.S / u.S)
        @. du.D_h = p.deb.k_D_h * (SL_max / SL) * (u.C_W - u.D_h) - u.D_h * (du.S / u.S)
    end
end

@inline function C_Wdot!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )
    du.C_W = zeros(length(u.C_W))
end

@inline function X_pdot!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )
    du.X_p = embryo(u.life_stage) ? 0. : p.glb.Xdot_in - du.I
end

@inline function X_embdot!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )
    du.X_emb = embryo(u.life_stage) ? -du.I : 0.0
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

    determine_life_stage!(du, u, p, t)

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