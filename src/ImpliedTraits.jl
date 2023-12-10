"""
Calculate maximum structural length slmax [m^(1/3)]
$(TYPEDSIGNATURES)
"""
function calc_slmax(anm::DEBParams)
    return ((anm.kappa * anm.Idotmax_rel * anm.eta_IA) / anm.k_M)
end

"""
Calculate maximum structural mass smax [m]
$(TYPEDSIGNATURES)
"""
function calc_smax(anm::DEBParams)
    return calc_slmax(anm)^3
end