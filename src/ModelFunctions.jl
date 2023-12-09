
"""
Determine the current life stage. 
$(TYPEDSIGNATURES)
"""
function determine_life_stage(
    H::Float64,
    H_p::Float64,
    X_emb::Float64,
    X_emb_int::Float64
)
    if H >= H_p
        return 5 # adult
    elseif H >= H_46
        return 4 # juvenile
    elseif H >= H_42
        return 3 # metamorph
    elseif isapprox(X_emb, 0, atol = RTOL_EMB * X_emb_int)
        return 2 # larva
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
function larva(life_stage::Int64)
    return life_stage == 2
end

"""
$(TYPEDSIGNATURES)
"""
function metamorph(life_stage::Int64)
    return life_stage == 3
end

"""
$(TYPEDSIGNATURES)
"""
function juvenile(life_stage::Int64)
    return life_stage == 4
end

"""
$(TYPEDSIGNATURES)
"""
function adult(life_stage::Int64)
    return life_stage == 5
end

function functional_response(X_V::Float64, K_X::Float64)
    return X_V / (X_V + K_X)
end

"""
Calculate ingestion rate for a single resource.
$(TYPEDSIGNATURES)
"""
function Idot(
    X_p::Float64,
    volume::Float64,
    X_emb::Float64,
    life_stage::Int64,
    S::Float64,
    Idot_max_rel::Float64,
    K_X::Float64 # TODO: replace this with K_X
    )
    
    if sum([juvenile(life_stage), adult(life_stage)]) > 0 # these life stages feed from external resource
        let X_p_V = X_p / volume, f_X = X_p_V / (X_p_V + K_X)
            return f_X * Idot_max_rel * S^(2/3) 
        end
    elseif embryo(life_stage) # embryos feed from the vitellus
        return min(X_emb, S^(2 / 3) * Idot_max_rel)
    else 
        error("Life stage not recognized: $(life_stage)")
    end
end


"""
Assimilation rate
$(TYPEDSIGNATURES)
"""
function Adot(Idot::Float64, eta_IA::Float64, Cdot::Float64, eta_SA::Float64)
    return Idot * eta_IA + Cdot * eta_SA
end


"""
Somatic maintenance rate
$(TYPEDSIGNATURES)
"""
function Mdot(S::Float64, k_M::Float64)
    return max(0, S * k_M)
end



"""
Maturity maintenance rate
$(TYPEDSIGNATURES)
"""
function Jdot(H::Float64, k_J::Float64)
    return H * k_J
end


"""
Positive somatic growth
$(TYPEDSIGNATURES)
"""
function Sdot_positive(eta_AS::Float64, kappa::Float64, Adot::Float64, Mdot::Float64)
    return eta_AS * (kappa * Adot - Mdot)
end

"""
Negative somatic growth
($(TYPEDSIGNATURES))
"""
function Sdot_negative(Mdot::Float64, eta_SA::Float64, kappa::Float64, Adot::Float64)
    return -(Mdot / eta_SA - kappa * Adot)
end

"""
Somatic growth rate
$(TYPEDSIGNATURES)
"""
function Sdot(
    Adot::Float64,
    Mdot::Float64,
    kappa::Float64,
    eta_AS::Float64,
    eta_SA::Float64,
    Cdot::Float64,
)
    let Sdot = (Sdot_positive(eta_AS, kappa, Adot, Mdot)) # calculate structural growth 
        if Sdot < 0 # if growth is negative, apply the shrinking equation
            Sdot = Sdot_negative(Mdot, eta_SA, kappa, Adot)
        end
        return Sdot - Cdot # return derivative with account for mobilized structure
    end
end

"""
Maturation rate. 
Maturity is dissipated energy and can therefore not be burned to cover maintenance costs. 
For simplicity, we currently assume that maturity maintenance will not be covered of the 1-kappa flux is insufficient.
$(TYPEDSIGNATURES)
"""
function Hdot(kappa::Float64, Adot::Float64, Jdot::Float64, adult::Bool)
    if adult == false
        return max(0, ((1 - kappa) * Adot) - Jdot)
    else
        return 0.0
    end
end

"""
Reproduction rate.
$(TYPEDSIGNATURES)
"""
function Rdot(kappa::Float64, Adot::Float64, Jdot::Float64, adult::Bool)
    if adult
        return (1 - kappa) * Adot - Jdot
    else
        return 0.0
    end
end


"""
Calculation of the total dissipation flux, equal to the sum of maintenance costs and overheads paid for assimilation, mobilization, growth and reproduction.
$(TYPEDSIGNATURES)
"""
function Qdot(
    Idot::Float64, 
    eta_IA::Float64, 
    kappa::Float64, 
    Sdot::Float64, 
    S::Float64,
    eta_AS::Float64, 
    eta_SA::Float64, 
    eta_AR::Float64,
    Mdot::Float64, 
    Jdot::Float64,
    Hdot::Float64
    )
    # dissipation fluxes for the individual processes
    let Qdot_A, Qdot_S, Qdot_C, Qdot_R
        Qdot_A = Idot * (1 - eta_IA)
        Qdot_S = Sdot >= 0 ? Sdot * (1 - eta_AS)/eta_AS : Sdot * (eta_SA - 1)
        Qdot_C = k_C * S * (1 - kappa) * (1 - eta_SA)
        Qdot_R = Rdot * (1 - eta_AR)/eta_AR
        return Qdot_A + Qdot_S + Qdot_C + Qdot_R + Mdot + Jdot + Hdot
    end
end
