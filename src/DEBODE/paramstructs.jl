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
`SpeciesParams` contain the species-level DEB and TKTD parameters. 
The default species parameters can be initialized using `SpeciesParams()`. 

The default parameters define an organism whose life-history is similar to 
*D. magna*.
"""
@with_kw mutable struct SpeciesParams <: AbstractSpeciesParams
    Z::Distribution = Dirac(1.) # agent variability is accounted for in the zoom factor. This can be set to a Dirac distribution if a zoom factor should be applied without introducing agent variability.
    propagate_zoom::NamedTuple = ( 
            Idot_max_rel_0 = true, 
            Idot_max_rel_emb_0 = true,
            X_emb_int_0 = true,
            H_p_0 = true, 
            K_X_0 = true
        )         
    T_A::Float64 = 8000. # Arrhenius temperature [K]
    T_ref::Float64 = 293.15 # reference temperature [K]
    X_emb_int_0::Float64 = 19.42 # initial vitellus [μgC]
    K_X_0::Float64 = 10_000 # half-saturation constant for food uptake [μgC L-1]
    Idot_max_rel_0::Float64 = 22.9 # maximum size-specific ingestion rate [μgC μgC^-(2/3) d-1]
    Idot_max_rel_emb_0::Float64 = 22.9 # size-specific embryonic ingestion rate [μgC μgC^-(2/3) d-1]
    kappa_0::Float64 = 0.539 # somatic allocation fraction [-]
    eta_IA_0::Float64 = 0.33 # assimilation efficiency [-]
    eta_AS_0::Float64 = 0.8 # growth efficiency [-]
    eta_SA::Float64 = 0.8 # shrinking efficiency [-]
    eta_AR_0::Float64 = 0.95 # reproduction efficiency [-]
    k_M_0::Float64 = 0.59 # somatic maintenance rate constant [d^-1]
    k_J_0::Float64 = 0.504 # maturity maintenance rate constant [d^-1]
    H_p_0::Float64 = 1/3 #100. # maturity at puberty [μgC]
    
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

    # the following are parameters which are currently only relevant for the  ABM
    
    f_Xthr::Float64 = 0.5 # functional response threshold for starvation mortality
    s_min::Float64 = 0.25 # daily survival mortality at complete food deprivation

    a_max = Truncated(Normal(60, 6), 0, Inf) # maximum life span 
    tau_R = 2. # reproduction interval
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

    agn.X_emb_int_0 = spc.X_emb_int_0 * agn.Z
    agn.H_p_0 = spc.H_p_0 * agn.Z
    agn.K_X_0 = spc.K_X_0 * agn.Z^(1/3) # K_X = Idot_max_rel / F_m with F_m = searching rate

    for param in fieldnames(typeof(spc.propagate_zoom)) # iterate over other parameters which may be affected by Z
        if getproperty(spc.propagate_zoom, param) # check whether propagation of Z should occur for this parameter
            setproperty!(agn, param, getproperty(spc, param) * agn.Z) # assign the agent value by adjusting the hyperparameter
        else # if Z should not be propagated to this parameter, 
            setproperty!(agn, param, getproperty(spc, param)) # set the agent-specific value equal to the population mean
        end
    end 
end

abstract type AbstractAgentParams <: AbstractParams end

"""
    AgentParams(spc::Union{AbstractParams,NamedTuple}) <: AbstractAgentParams

DEBBase defines parameters on three levels: Global, species-level and agent-level. 

Users typically do not need to interact with `AgentParams`. 
However, if you are building a modification of the base model and you wish to change 
how individual variability is treated, 
you need to define your own `AbstractAgentParams`. 
Note that `DEBODE.simulator` also accepts an `AgentParamType` as keyword argument.

AgentParams holds the parameters which are subject to agent variability. 
This is in contrast to SpeciesParams, which define parameters on the species-level.


Implementing a custom model might require to modify `AgentParams` if individual 
    variability should occur in other variables.
"""
@with_kw mutable struct AgentParams <: AbstractAgentParams
    Z::Float64
    Idot_max_rel_0::Float64
    Idot_max_rel_emb_0::Float64
    X_emb_int_0::Float64
    H_p_0::Float64
    K_X_0::Float64
    a_max::Float64
    
    """
    Initialize AgentParams from SpeciesParams `spc`.
    """
    function AgentParams(spc::Union{AbstractParams,NamedTuple})
        agn = new()
        individual_variability!(agn, spc)
        return agn
    end
end

"""
A `Params` object contains global parameters `glb` and species parameters `spc` 
(including TKTD-parameters). \\
Initialize the default parameter collection with `Params()`.
"""
@with_kw mutable struct Params <: AbstractParamCollection
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


