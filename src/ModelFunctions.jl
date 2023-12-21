
"""
Determine the current life stage. 
$(TYPEDSIGNATURES)
"""
function determine_life_stage(
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
function embryo(life_stage::Int64)
    return life_stage == 1
end


"""
$(TYPEDSIGNATURES)
"""
function juvenile(life_stage::Int64)
    return life_stage == 2
end

"""
$(TYPEDSIGNATURES)
"""
function adult(life_stage::Int64)
    return life_stage == 3
end


function functional_response(deb::AbstractParams; X_V::Float64)
    return X_V / (X_V + deb.K_X)
end

"""
Calculate ingestion rate for a single resource.
$(TYPEDSIGNATURES)
"""
function Idot(
    glb::AbstractParams,
    deb::AbstractParams;
    X_p::Float64,
    life_stage::Int64,
    S::Float64
    )
    
    if life_stage > 1 # juveniles and adults feed from external resource
        let X_p_V = X_p / glb.V_patch, 
            f_X = functional_response(deb; X_V = X_p_V)
            return f_X * deb.Idot_max_rel_emb * S^(2/3) 
        end
    else # embryos feed from the vitellus
        return S^(2 / 3) * deb.Idot_max_rel
    end
end

"""
Assimilation rate
$(TYPEDSIGNATURES)
"""
function Adot(deb::AbstractParams; Idot::Float64)
    return Idot * deb.eta_IA
end

"""
Somatic maintenance rate
$(TYPEDSIGNATURES)
"""
function Mdot(deb::AbstractParams; S::Float64)
    return max(0, S * deb.k_M)
end

"""
Maturity maintenance rate
$(TYPEDSIGNATURES)
"""
function Jdot(deb::AbstractParams; H::Float64)
    return H * deb.k_J
end

"""
Positive somatic growth
$(TYPEDSIGNATURES)
"""
function Sdot_positive(deb::AbstractParams; Adot::Float64, Mdot::Float64)
    return deb.eta_AS * (deb.kappa * Adot - Mdot)
end

"""
Negative somatic growth
($(TYPEDSIGNATURES))
"""
function Sdot_negative(deb::AbstractParams; Mdot::Float64, Adot::Float64)
    return -(Mdot / deb.eta_SA - deb.kappa * Adot)
end

"""
Somatic growth rate
$(TYPEDSIGNATURES)
"""
function Sdot(
    deb::AbstractParams;
    Mdot::Float64,
    Adot::Float64
)
    let Sdot = (Sdot_positive(deb; Mdot = Mdot, Adot = Adot)) # calculate structural growth 
        if Sdot < 0 # if growth is negative, apply the shrinking equation
            Sdot = Sdot_negative(deb; Mdot = Mdot, Adot = Adot)
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
function Hdot(deb::AbstractParams; Adot::Float64, Jdot::Float64, adult::Bool)
    if adult == false
        return max(0, ((1 - deb.kappa) * Adot) - Jdot)
    else
        return 0.0
    end
end

"""
Reproduction rate.
$(TYPEDSIGNATURES)
"""
function Rdot(deb::AbstractParams; Adot::Float64, Jdot::Float64, adult::Bool)
    if adult
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
Definition of reserveless DEB derivatives. 
$(TYPEDSIGNATURES)
"""
function DEB!(du, u, p, t)
    glb, deb = p
    # unpack state variables
    X_p, X_emb, S, H, R = u

    S = max(0, S) # control for negative values

    life_stage = determine_life_stage(deb; H = H, X_emb = X_emb)
    
    idot = Idot(glb, deb; X_p = X_p, life_stage = life_stage, S = S)
    xdot = embryo(life_stage) ? 0. : glb.Xdot_in - idot
    xembdot = embryo(life_stage) ? -idot : 0.0
    adot = Adot(deb; Idot = idot, )
    mdot = Mdot(deb; S = S)
    jdot = Jdot(deb; H = H)
    sdot = Sdot(deb; Adot = adot, Mdot = mdot) 
    hdot = Hdot(deb; Adot = adot, Jdot = jdot, adult = adult(life_stage))
    rdot = Rdot(deb; Adot = adot, Jdot = jdot, adult = adult(life_stage))

    # update du/dt
    du[01] = xdot
    du[02] = xembdot
    du[03] = sdot
    du[04] = hdot
    du[05] = rdot
end


"""
TK for DEB-TKTD model, including effect of surface area to volume ratio and dilution by growth.
$(TYPEDSIGNATURES)
"""
function DEBTK(
    tktd::TKTDParams; 
    C_W::Float64, 
    D::Float64, 
    S::Float64, 
    Sdot::Float64, 
    S_max::Float64
    )
    let AV_ratio = S^(2/3) / S, AV_max_ratio = S_max^(2/3) / S_max
        return tktd.k_D * (AV_max_ratio / AV_ratio) * (C_W - D) - D * (Sdot/S)
    end
end

