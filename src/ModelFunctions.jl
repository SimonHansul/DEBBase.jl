
"""
Determine the current life stage. 
$(TYPEDSIGNATURES)
"""
function determine_life_stage(
    p::DEBParams;
    H::Float64,
    X_emb::Float64,
)
    if H >= p.H_p
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


function functional_response(p::DEBParams; X_V::Float64)
    return X_V / (X_V + p.K_X)
end

"""
Calculate ingestion rate for a single resource.
$(TYPEDSIGNATURES)
"""
function Idot(
    g::GlobalParams,
    p::DEBParams;
    X_p::Float64,
    X_emb::Float64,
    life_stage::Int64,
    S::Float64
    )
    
    if sum([juvenile(life_stage), adult(life_stage)]) > 0 # these life stages feed from external resource
        let X_p_V = X_p / g.volume, f_X = X_p_V / (X_p_V + p.K_X)
            return f_X * p.Idot_max_rel * S^(2/3) 
        end
    elseif embryo(life_stage) # embryos feed from the vitellus
        return min(X_emb, S^(2 / 3) * p.Idot_max_rel)
    else 
        error("Life stage not recognized: $(life_stage)")
    end
end

"""
Assimilation rate
$(TYPEDSIGNATURES)
"""
function Adot(p::DEBParams; Idot::Float64)
    return Idot * p.eta_IA
end

"""
Somatic maintenance rate
$(TYPEDSIGNATURES)
"""
function Mdot(p::DEBParams; S::Float64)
    return max(0, S * p.k_M)
end

"""
Maturity maintenance rate
$(TYPEDSIGNATURES)
"""
function Jdot(p::DEBParams; H::Float64)
    return H * p.k_J
end

"""
Positive somatic growth
$(TYPEDSIGNATURES)
"""
function Sdot_positive(p::DEBParams; Adot::Float64, Mdot::Float64)
    return p.eta_AS * (p.kappa * Adot - Mdot)
end

"""
Negative somatic growth
($(TYPEDSIGNATURES))
"""
function Sdot_negative(p::DEBParams; Mdot::Float64, Adot::Float64)
    return -(Mdot / p.eta_SA - p.kappa * Adot)
end

"""
Somatic growth rate
$(TYPEDSIGNATURES)
"""
function Sdot(
    p::DEBParams;
    Mdot::Float64,
    Adot::Float64
)
    let Sdot = (Sdot_positive(p; Mdot = Mdot, Adot = Adot)) # calculate structural growth 
        if Sdot < 0 # if growth is negative, apply the shrinking equation
            Sdot = Sdot_negative(p; Mdot = Mdot, Adot = Adot)
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
function Hdot(p::DEBParams; Adot::Float64, Jdot::Float64, adult::Bool)
    if adult == false
        return max(0, ((1 - p.kappa) * Adot) - Jdot)
    else
        return 0.0
    end
end

"""
Reproduction rate.
$(TYPEDSIGNATURES)
"""
function Rdot(p::DEBParams; Adot::Float64, Jdot::Float64, adult::Bool)
    if adult
        return (1 - p.kappa) * Adot - Jdot
    else
        return 0.0
    end
end


"""
Calculation of the total dissipation flux, equal to the sum of maintenance costs and overheads paid for assimilation, mobilization, growth and reproduction.
$(TYPEDSIGNATURES)
"""
function Qdot(
    p::DEBParams;
    Idot::Float64, 
    Sdot::Float64, 
    Mdot::Float64, 
    Jdot::Float64,
    Hdot::Float64
    )
    # dissipation fluxes for the individual processes
    let Qdot_A, Qdot_S, Qdot_C, Qdot_R
        Qdot_A = Idot * (1 - p.eta_IA)
        Qdot_S = Sdot >= 0 ? Sdot * (1 - p.eta_AS)/p.eta_AS : Sdot * (p.eta_SA - 1)
        Qdot_R = Rdot * (1 - p.eta_AR)/p.eta_AR
        return Qdot_A + Qdot_S + Qdot_C + Qdot_R + Mdot + Jdot + Hdot
    end
end


"""
Definition of reserveless DEB derivatives. 
$(TYPEDSIGNATURES)
"""
function DEB!(du, u, p, t)
    glb, anm = p
    # unpack state variables
    X_p, X_emb, S, H, R = u

    S = max(0, S) # control for negative values

    life_stage = determine_life_stage(anm; H = H, X_emb = X_emb)
   
    idot = Idot(glb, anm; X_p = X_p, X_emb = X_emb, life_stage = life_stage, S = S)
    xdot = embryo(life_stage) ? 0. : glb.Xdot_in - idot
    xembdot = embryo(life_stage) ? -idot : 0.0
    adot = Adot(anm; Idot = idot, )
    mdot = Mdot(anm; S = S)
    jdot = Jdot(anm; H = H)
    sdot = Sdot(anm; Adot = adot, Mdot = mdot) 
    hdot = Hdot(anm; Adot = adot, Jdot = jdot, adult = adult(life_stage))
    rdot = Rdot(anm; Adot = adot, Jdot = jdot, adult = adult(life_stage))

    # update du/dt
    du[01] = xdot
    du[02] = xembdot
    du[03] = sdot
    du[04] = hdot
    du[05] = rdot
end
