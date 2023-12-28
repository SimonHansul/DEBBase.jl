
"""
Determine the current life stage. 
$(TYPEDSIGNATURES)
"""
@inline function determine_life_stage(
    deb::AbstractParams;
    H::Float64,
    X_emb::Float64,
)
    if H >= deb.H_p
        return 3 # adult
    elseif X_emb <= 0
        return 2 # juvenile
    else
        return 1 # embryo
    end
end


"""
$(TYPEDSIGNATURES)
"""
@inline function embryo(life_stage::Int64)
    return life_stage == 1
end


"""
$(TYPEDSIGNATURES)
"""
@inline function juvenile(life_stage::Int64)
    return life_stage == 2
end

"""
$(TYPEDSIGNATURES)
"""
@inline function adult(life_stage::Int64)
    return life_stage == 3
end


@inline function functional_response(deb::AbstractParams; X_V::Float64)
    return X_V / (X_V + deb.K_X)
end

"""
Calculate ingestion rate for a single resource.
$(TYPEDSIGNATURES)
"""
@inline function Idot(
    glb::AbstractParams,
    deb::AbstractParams,
    du::ComponentArray,
    u::ComponentArray;
    life_stage::Int64
    )
    
    if life_stage > 1 # juveniles and adults feed from external resource
        let X_p_V = u.X_p / glb.V_patch, 
            f_X = functional_response(deb; X_V = X_p_V)
            return f_X * deb.Idot_max_rel_emb * u.S^(2/3) 
        end
    else # embryos feed from the vitellus
        return u.S^(2 / 3) * deb.Idot_max_rel
    end
end

"""
Assimilation rate
$(TYPEDSIGNATURES)
"""
@inline function Adot(
    glb::AbstractParams,
    deb::AbstractParams,
    du::ComponentArray,
    u::ComponentArray;
    Idot::Float64
    )
    return Idot * deb.eta_IA
end

"""
Somatic maintenance rate
$(TYPEDSIGNATURES)
"""
@inline function Mdot(
    glb::AbstractParams,
    deb::AbstractParams,
    du::ComponentArray,
    u::ComponentArray;)
    return max(0, u.S * deb.k_M)
end

"""
Maturity maintenance rate
$(TYPEDSIGNATURES)
"""
@inline function Jdot(
    glb::AbstractParams,
    deb::AbstractParams,
    du::ComponentArray,
    u::ComponentArray
    )
    return u.H * deb.k_J
end

"""
Positive somatic growth
$(TYPEDSIGNATURES)
"""
@inline function Sdot_positive(
    glb::AbstractParams,
    deb::AbstractParams,
    du::ComponentArray,
    u::ComponentArray;
    Adot::Float64,
    Mdot::Float64
    )
    return deb.eta_AS * (deb.kappa * Adot - Mdot)
end

"""
Negative somatic growth
($(TYPEDSIGNATURES))
"""
@inline function Sdot_negative(
    glb::AbstractParams,
    deb::AbstractParams,
    du::ComponentArray,
    u::ComponentArray;
    Adot::Float64,
    Mdot::Float64
    )
    return -(Mdot / deb.eta_SA - deb.kappa * Adot)
end

"""
Somatic growth rate
$(TYPEDSIGNATURES)
"""
@inline function Sdot(
    glb::AbstractParams,
    deb::AbstractParams,
    du::ComponentArray,
    u::ComponentArray;
    Adot::Float64,
    Mdot::Float64
    )
    let Sdot = (Sdot_positive(glb, deb, du, u; Mdot = Mdot, Adot = Adot)) # calculate structural growth 
        if Sdot < 0 # if growth is negative, apply the shrinking equation
            Sdot = Sdot_negative(glb, deb, du, u; Mdot = Mdot, Adot = Adot)
        end
        return Sdot # return derivative with account for mobilized structure
    end
end

"""
Maturation rate. 
Maturity is dissipated energy and can therefore not be burned to cover maintenance costs. 
For simplicity, we currently assume that maturity maintenance will not be covered of the 1-kappa flux is insufficient.
$(TYPEDSIGNATURES)
"""
@inline function Hdot(
    glb::AbstractParams,
    deb::AbstractParams,
    du::ComponentArray,
    u::ComponentArray; 
    Adot::Float64, 
    Jdot::Float64, 
    life_stage::Int64
    )
    if adult(life_stage) == false
        return max(0, ((1 - deb.kappa) * Adot) - Jdot)
    else
        return 0.0
    end
end

"""
Reproduction rate.
$(TYPEDSIGNATURES)
"""
@inline function Rdot(
    glb::AbstractParams,
    deb::AbstractParams,
    du::ComponentArray,
    u::ComponentArray; 
    Adot::Float64, 
    Jdot::Float64,
    life_stage::Int64
    )
    if adult(life_stage)
        return (1 - deb.kappa) * Adot - Jdot
    else
        return 0.0
    end
end

"""
Calculation of the total dissipation flux, equal to the sum of maintenance costs and overheads paid for assimilation, mobilization, growth and reproduction.
$(TYPEDSIGNATURES)
"""
function Qdot(
    deb::AbstractParams;
    Idot::Float64, 
    Sdot::Float64, 
    Mdot::Float64, 
    Jdot::Float64,
    Hdot::Float64
    )
    # dissipation fluxes for the individual processes
    let Qdot_A, Qdot_S, Qdot_C, Qdot_R
        Qdot_A = Idot * (1 - deb.eta_IA)
        Qdot_S = Sdot >= 0 ? Sdot * (1 - deb.eta_AS)/deb.eta_AS : Sdot * (deb.eta_SA - 1)
        Qdot_R = Rdot * (1 - deb.eta_AR)/deb.eta_AR
        return Qdot_A + Qdot_S + Qdot_C + Qdot_R + Mdot + Jdot + Hdot
    end
end

"""
TK for DEB-TKTD model, including effect of surface area to volume ratio and dilution by growth. 
If `D` and is given as a Vector, TK is only stressor-specific but not PMoA-specific. 
$(TYPEDSIGNATURES)
"""
@inline function Ddot(
    glb::AbstractParams,
    deb::AbstractParams,
    du::ComponentVector,
    u::ComponentVector
    )

    let L_S = u.S^(1/3) / u.S, # strucutral length (g^(1/3))
        # TODO: move calculation of L_S_max out so it is only calculated once, not at every step
        L_S_max = calc_S_max(deb)^(1/3) # maximum structural length (g^(1/3))
        return @. deb.k_D * (L_S_max / L_S) * (u.C_W - u.D) - u.D * du.S / u.S
    end
end



"""
Definition of reserveless DEB derivatives. 
$(TYPEDSIGNATURES)
"""
function DEB!(du, u, p, t)

    #### boilerplate

    glb::GlobalBaseParams, deb::DEBBaseParams = p # unpack parameters
    @unpack X_p, X_emb, S, H, R, D, C_W = u # unpack state variables

    S = max(0, S) # control for negative values
    life_stage = determine_life_stage(deb; H = H, X_emb = X_emb)
    
    #### auxiliary state variables

    idot = Idot(glb, deb, du, u; life_stage = life_stage)
    adot = Adot(glb, deb, du, u; Idot = idot)
    mdot = Mdot(glb, deb, du, u)
    jdot = Jdot(glb, deb, du, u)

    #### major state variables

    du.X_p = embryo(life_stage) ? 0. : glb.Xdot_in - idot
    du.X_emb = embryo(life_stage) ? -idot : 0.0
    du.S = Sdot(glb, deb, du, u; Adot = adot, Mdot = mdot)
    du.H = Hdot(glb, deb, du, u; Adot = adot, Jdot = jdot, life_stage = life_stage)
    du.R = Rdot(glb, deb, du, u; Adot = adot, Jdot = jdot, life_stage = life_stage)
    du.D = Ddot(glb, deb, du, u)
    du.C_W = zeros(length(glb.C_W))
end