"""
`GlobalParams` contain the global parameters (simulated timespan `t_max`, nutrient influx rate `Xdot_in`, etc.)
"""
@with_kw mutable struct GlobalParams <: AbstractGlobalParams
    N0::Int64 = 1 #  initial number of individuals [#]
    t_max::Float64 = 21. # maximum simulation time [d]
    Xdot_in::Float64 = 1200. # nutrient influx set to be a little above the absolute maximum ingestion rate according to default SpeciesParams [μg C d-1]
    k_V::Float64 = 0. # chemostat dilution rate [d-1]
    V_patch::Float64 = 0.05 # volume of a patch (L) (or the entire similated environment) [L]
    C_W::Vector{Float64} = [0.] # external chemical concentrations [μg L-1]
    T::Float64 = 293.15 # ambient temperature [K]
end

"""
`SpeciesParams` contain population means of DEB and TKTD parameters. 
Default values are roughly reproduce life-history of Daphnia minda (model currency mass carbon in μgC) exposued to Azoxystrobin (mg/L). \\
DEBBase.jl uses a hierarchical modelling approach where the `SpeciesParams` are  parameters which are common across all agents of a species, 
and `ODEAgentParams` contain parameters which are specific for an individual. \\
Variability is given by the zoom factor `Z::Distribution`, which is always applied to the surface-area specific ingestion rates 
and can optionally propagate to parameters indicated in `propagate_zoom::NTuple`. \\
`Z` is `Dirac(1)` by default, i.e. there is no agent variability in the default parameters. \\
"""
@with_kw mutable struct SpeciesParams <: AbstractSpeciesParams
    Z::Distribution = Dirac(1.) # agent variability is accounted for in the zoom factor. This can be set to a Dirac distribution if a zoom factor should be applied without introducing agent variability.
    propagate_zoom::@NamedTuple{X_emb_int_0::Bool, H_p_0::Bool, K_X_0::Bool} = (
        X_emb_int_0 = true,
        H_p_0 = true, 
        K_X_0 = true
        ) # Parameters to which Z will be propagated. Z is *always* applied to `Idot_max_rel_0` (with appropriate scaling).
        
    T_A::Float64 = 8000. # Arrhenius temperature [K]
    T_ref::Float64 = 293.15 # reference temperature [K]
    X_emb_int_0::Float64 = 19.42 # initial vitellus [μgC]
    K_X_0::Float64 = 1. # half-saturation constant for food uptake [μgC L-1]
    Idot_max_rel_0::Float64 = 22.9 # maximum size-specific ingestion rate [μgC μgC^-(2/3) d-1]
    Idot_max_rel_emb_0::Float64 = 22.9 # size-specific embryonic ingestion rate [μgC μgC^-(2/3) d-1]
    kappa_0::Float64 = 0.539 # somatic allocation fraction [-]
    eta_IA_0::Float64 = 0.33 # assimilation efficiency [-]
    eta_AS_0::Float64 = 0.8 # growth efficiency [-]
    eta_SA::Float64 = 0.8 # shrinking efficiency [-]
    eta_AR_0::Float64 = 0.95 # reproduction efficiency [-]
    k_M_0::Float64 = 0.59 # somatic maintenance rate constant [d^-1]
    k_J_0::Float64 = 0.504 # maturity maintenance rate constant [d^-1]
    H_p_0::Float64 = 100. # maturity at puberty [μgC]
    
    k_D_G::Vector{Float64} = [0.00] # toxicokinetic rate constants | PMoA growth efficiency
    k_D_M::Vector{Float64} = [0.00] # toxicokinetic rate constants | PMoA maintenance costs
    k_D_A::Vector{Float64} = [0.00] # toxicokinetic rate constants | PMoA assimilation efficiency
    k_D_R::Vector{Float64} = [0.38] # toxicokinetic rate constants | PMoA reproduction efficiency
    k_D_h::Vector{Float64} = [0.00] # toxicokinetic rate constants | PMoA hazard rate
    
    drc_functs_G::Vector{Function} = [LL2]  # dose-response functions | PMoA growth efficiency
    drc_functs_M::Vector{Function} = [LL2M] # dose-response functions | PMoA maintenance costs
    drc_functs_A::Vector{Function} = [LL2]  # dose-response functions | PMoA assimilation efficiency
    drc_functs_R::Vector{Function} = [LL2]  # dose-response functions | PMoA reproduction efficiency
    drc_functs_h::Vector{Function} = [LL2h] # dose-response functions | PMoA hazard rate

    e_G::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA growth efficiency
    e_M::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA maintenance costs
    e_A::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA assimilation efficiency
    e_R::Union{Nothing,Vector{Float64}} = [167.] # sensitivity parameters | PMoA reproduction efficiency
    e_h::Union{Nothing,Vector{Float64}} = [1e10] # sensitivity parameters | PMoA hazard rate

    b_G::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA growth efficiency
    b_M::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA maintenance costs
    b_A::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA assimilation efficiency
    b_R::Union{Nothing,Vector{Float64}} = [0.93] # slope parameters | PMoA reproduction efficiency
    b_h::Union{Nothing,Vector{Float64}} = [1e10] # slope parameters | PMoA reproduction efficiency
    
    a_max = Truncated(Normal(60, 6), 0, Inf) # maximum life span - currently only used by ABM
    tau_R = 2. # reproduction interval - currently only used by ABM
end

"""
    individual_variability!(p::Ref{AbstractParams})

Induce agent variability in spc parameters via zoom factor `Z`. 
`Z` is sampled from the corresponding distribution given in `p` and assumed to represent a ratio between maximum structurel *masses* (not lengths), 
so that the surface area-specific ingestion rate `Idot_max_rel_0` scales with `Z^(1/3)` and parameters which represent masses or energy pools scales with `Z`.
"""
function individual_variability!(agn::Union{AbstractParams,NamedTuple}, spc::Union{AbstractParams,NamedTuple})
    agn.Z = rand(spc.Z) # sample zoom factor Z for agent from distribution
    agn.a_max = rand(spc.a_max) # sample maximum age from distribution - currently only used by ABM

    agn.Idot_max_rel_0 = spc.Idot_max_rel_0 * agn.Z^(1/3) # Z is always applied to Idot_max_rel_0
    agn.Idot_max_rel_emb_0 = spc.Idot_max_rel_emb_0 * agn.Z^(1/3) #, including the value for embryos

    for param in fieldnames(typeof(spc.propagate_zoom)) # iterate over other parameters which may be affected by Z
        if getproperty(spc.propagate_zoom, param) # check whether propagation of Z should occur for this parameter
            setproperty!(agn, param, getproperty(spc, param) * agn.Z) # assign the agent value by adjusting the hyperparameter
        else # if Z should not be propagated to this parameter, 
            setproperty!(agn, param, getproperty(spc, param)) # set the agent-specific value equal to the population mean
        end
    end 
end

"""
    ODEAgentParams(spc::Union{AbstractParams,NamedTuple})
ODEAgentParams are subject to agent variability. 
This is in contrast to SpeciesParams, which define parameters on the species-level.
"""
@with_kw mutable struct ODEAgentParams <: AbstractParams
    Z::Float64
    Idot_max_rel_0::Float64
    Idot_max_rel_emb_0::Float64
    X_emb_int_0::Float64
    H_p_0::Float64
    K_X_0::Float64
    
    """
    Initialize ODEAgentParams from SpeciesParams `spc`.
    """
    function ODEAgentParams(spc::Union{AbstractParams,NamedTuple})
        agn = new()
        individual_variability!(agn, spc)
        return agn
    end
end

"""
A `DEBParamCollection` contains global parameters `glb` and spc parameters `spc` (including TKTD-parameters). \\
Initialize the default parameter collection with `DEBParamCollection()`.
"""
@with_kw mutable struct DEBParamCollection <: AbstractParamCollection
    glb::Union{AbstractParams,NamedTuple} = GlobalParams()
    spc::Union{AbstractParams,NamedTuple} = SpeciesParams()
    agn::Union{Nothing,AbstractParams} = nothing

    @assert length(spc.k_D_G) >= length(glb.C_W) "Length of k_D_G is not at least length of C_W"
    @assert length(spc.k_D_M) >= length(glb.C_W) "Length of k_D_M is not at least length of C_W"
    @assert length(spc.k_D_A) >= length(glb.C_W) "Length of k_D_A is not at least length of C_W"
    @assert length(spc.k_D_R) >= length(glb.C_W) "Length of k_D_R is not at least length of C_W"
    @assert length(spc.k_D_h) >= length(glb.C_W) "Length of k_D_h is not at least length of C_W"
    @assert length(spc.drc_functs_G) >= length(glb.C_W) "Length of drc_functs_G is not at least length of C_W"
    @assert length(spc.drc_functs_M) >= length(glb.C_W) "Length of drc_functs_M is not at least length of C_W"
    @assert length(spc.drc_functs_A) >= length(glb.C_W) "Length of drc_functs_A is not at least length of C_W"
    @assert length(spc.drc_functs_R) >= length(glb.C_W) "Length of drc_functs_R is not at least length of C_W"
    @assert length(spc.drc_functs_h) >= length(glb.C_W) "Length of drc_functs_h is not at least length of C_W"
    @assert length(spc.e_G) >= length(glb.C_W) "Length of e_G is not at least length of C_W"
    @assert length(spc.e_M) >= length(glb.C_W) "Length of e_M is not at least length of C_W"
    @assert length(spc.e_A) >= length(glb.C_W) "Length of e_A is not at least length of C_W"
    @assert length(spc.e_R) >= length(glb.C_W) "Length of e_R is not at least length of C_W"
    @assert length(spc.e_h) >= length(glb.C_W) "Length of e_h is not at least length of C_W"
    @assert length(spc.b_G) >= length(glb.C_W) "Length of b_G is not at least length of C_W"
    @assert length(spc.b_M) >= length(glb.C_W) "Length of b_M is not at least length of C_W"
    @assert length(spc.b_A) >= length(glb.C_W) "Length of b_A is not at least length of C_W"
    @assert length(spc.b_R) >= length(glb.C_W) "Length of b_R is not at least length of C_W"
    @assert length(spc.b_h) >= length(glb.C_W) "Length of b_h is not at least length of C_W"
end


