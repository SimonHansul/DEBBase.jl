"""
Calculate maximum structural length slmax [m^(1/3)]
$(TYPEDSIGNATURES)
"""
function calc_SL_max(spc::AbstractParams)::Float64
    return ((spc.kappa * spc.Idot_max_rel * spc.eta_IA) / spc.k_M)
end

"""
Calculate maximum structural mass smax [m]
$(TYPEDSIGNATURES)
"""
function calc_S_max(spc::AbstractParams)::Float64
    return calc_SL_max(spc)^3
end

"""
Set the maturity maintenance rate constant, 
assuming that the cumulative investment into maturity maintenance 
equals the cumulative investment into somatic maintenance.
$(TYPEDSIGNATURES)
"""
function k_J!(spc::AbstractParams)::nothing
    spc.k_J = ((1 - spc.kappa) / spc.kappa) * spc.k_M
end