"""
Calculate maximum structural length slmax [m^(1/3)]
$(TYPEDSIGNATURES)
"""
function calc_SL_max(deb::AbstractParams)
    return ((deb.kappa * deb.Idot_max_rel * deb.eta_IA) / deb.k_M)
end

"""
Calculate maximum structural mass smax [m]
$(TYPEDSIGNATURES)
"""
function calc_S_max(deb::AbstractParams)
    return calc_SL_max(deb)^3
end