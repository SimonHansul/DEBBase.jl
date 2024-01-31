"""
Calculate maximum structural length slmax [m^(1/3)]
$(TYPEDSIGNATURES)
"""
function calc_SL_max(deb::AbstractParams)::Float64
    return ((deb.kappa * deb.Idot_max_rel * deb.eta_IA) / deb.k_M)
end

"""
Calculate maximum structural mass smax [m]
$(TYPEDSIGNATURES)
"""
function calc_S_max(deb::AbstractParams)::Float64
    return calc_SL_max(deb)^3
end

"""
Set the maturity maintenance rate constant, 
assuming that the cumulative investment into maturity maintenance 
equals the cumulative investment into somatic maintenance.
$(TYPEDSIGNATURES)
"""
function k_J!(deb::AbstractParams)::nothing
    deb.k_J = ((1 - deb.kappa) / deb.kappa) * deb.k_M
end