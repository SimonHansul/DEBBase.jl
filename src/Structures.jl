abstract type AbstractSpeciesParams <: AbstractParams end
abstract type AbstractAgentParams <: AbstractParams end

#=
## General structures
=#

"""
`GlobalParams` contain the global parameters (simulated timespan `t_max`, nutrient influx rate `Xdot_in`, etc.)
"""
@with_kw mutable struct GlobalParams <: AbstractParams
    N0::Int64 = 1 #  initial number of individuals
    t_max::Float64 = 21. # maximum simulation time (d)
    Xdot_in::Float64 = 1200. # resource inflow rate
    k_V::Float64 = 0.1 # chemostat dilution rate
    V_patch::Float64 = 0.05 # volume of a patch (L) (or the entire similated environment)
    C_W::Vector{Float64} = [0.] # external chemical concentrations
    AgentType::DataType = DEBAgent # the type of agent simulated
    recordagentvars::Bool = true # record agent-level output?
    saveat::Float64 = 1. # when to save output (d)
end

"""
`SpeciesParams` contain population means of DEB and TKTD parameters. Default values are for Daphnia magna and Azoxystrobin in μg C. <br>
DEBBase.jl uses a hierarchical modelling approach where the `SpeciesParams` are  parameters which are common across all agents of a species, 
and `AgentParams` contain parameters which are specific for a species. <br>
Variability is given by the zoom factor `Z::Distribution`, which is always applied to the surface-area specific ingestion rates 
and can optionally propagate to parameters indicated in `propagate_zoom::NTuple`. <br>
`Z` is `Dirac(1)` by default, i.e. there is no agent variability in the default parameters. <br>
"""
@with_kw mutable struct SpeciesParams <: AbstractSpeciesParams
    Z::Distribution = Dirac(1.) # agent variability is accounted for in the zoom factor. This can be set to a Dirac distribution if a zoom factor should be applied without introducing agent variability.
    propagate_zoom::@NamedTuple{X_emb_int::Bool, H_p::Bool, K_X::Bool} = (X_emb_int = true, H_p = true, K_X = true) # Parameters to which Z will be propagated. Z is *always* applied to `Idot_max_rel` (with appropriate scaling).
    X_emb_int::Float64 = 19.42 # initial vitellus
    K_X::Float64 = 1. # half-saturation constant for food uptake
    Idot_max_rel::Float64 = 22.9 # maximum size-specific ingestion rate
    Idot_max_rel_emb::Float64 = 22.9 # size-specific embryonic ingestion rate
    kappa::Float64 = 0.539 # Somatic allocation fraction
    eta_IA::Float64 = 0.33 # Assimilation efficiency
    eta_AS::Float64 = 0.8 # Growth efficiency
    eta_SA::Float64 = 0.8 # Shrinking efficiency
    eta_AR::Float64 = 0.95 # Reproduction efficiency
    k_M::Float64 = 0.59 # Somatic maintenance rate constant
    k_J::Float64 = 0.504 # Maturity maintenance rate constant
    H_p::Float64 = 100. # Agent maturity at puberty
    
    k_D_G::Vector{Float64} = [0.] # Dominant rate constants | PMoA growth efficiency
    k_D_M::Vector{Float64} = [0.] # Dominant rate constants | PMoA maintenance costs
    k_D_A::Vector{Float64} = [0.] # Dominant rate constants | PMoA assimilation efficiency
    k_D_R::Vector{Float64} = [0.38] # Dominant rate constants | PMoA reproduction efficiency
    k_D_h::Vector{Float64} = [0.] # Dominant rate constants | PMoA hazard rate
    
    drc_functs_G::Vector{Function} = [LL2] # Dose-response functions | PMoA growth efficiency
    drc_functs_M::Vector{Function} = [LL2M] # Dose-response functions | PMoA maintenance costs
    drc_functs_A::Vector{Function} = [LL2] # Dose-response functions | PMoA assimilation efficiency
    drc_functs_R::Vector{Function} = [LL2] # Dose-response functions | PMoA reproduction efficiency
    drc_functs_h::Vector{Function} = [LL2h] # Dose-response functions | PMoA hazard rate

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

    c_G::Union{Nothing,Vector{Float64}} = nothing # placeholders for additional TD parameters. only relevant if custom drc_functs with more than two parameters are defined
    c_M::Union{Nothing,Vector{Float64}} = nothing
    c_A::Union{Nothing,Vector{Float64}} = nothing
    c_R::Union{Nothing,Vector{Float64}} = nothing
    c_h::Union{Nothing,Vector{Float64}} = nothing

    d_G::Union{Nothing,Vector{Function}} = nothing
    d_M::Union{Nothing,Vector{Function}} = nothing
    d_A::Union{Nothing,Vector{Function}} = nothing
    d_R::Union{Nothing,Vector{Function}} = nothing
    d_h::Union{Nothing,Vector{Function}} = nothing

    aux::Any = nothing # container for auxiliaray pramaters -- good for development purposes
end

"""
    agent_variability!(p::Ref{AbstractParams})
Induce agent variability in spc parameters via zoom factor `Z`. 
`Z` is sampled from the corresponding distribution given in `p` and assumed to represent a ratio between maximum structurel *masses* (not lengths), 
so that the surface area-specific ingestion rate `Idot_max_rel` scales with `Z^(1/3)` and parameters which represent masses or energy pools scales with `Z`.
"""
function agent_variability!(agn::AGN, spc::SPC) where {AGN <: AbstractParams, SPC <: AbstractParams}
    agn.Z = rand(spc.Z) # sample zoom factor Z for agent from distribution
    agn.Idot_max_rel = spc.Idot_max_rel * agn.Z^(1/3) # Z is always applied to Idot_max_rel
    agn.Idot_max_rel_emb = spc.Idot_max_rel_emb * agn.Z^(1/3) #, including the value for embryos

    for param in fieldnames(typeof(spc.propagate_zoom)) # iterate over other parameters which may be affected by Z
        if getproperty(spc.propagate_zoom, param) # check whether propagation of Z should occur for this parameter
            setproperty!(agn, param, getproperty(spc, param) * agn.Z) # assign the agent value by adjusting the hyperparameter
        else # if Z should not be propagated to this parameter, 
            setproperty!(agn, param, getproperty(spc, param)) # set the agent-specific value equal to the population mean
        end
    end 
end

"""
    AgentParams(spc::AbstractParams)
AgentParams are subject to agent variability. 
This is in contrast to SpeciesParams, which define parameters on the species-level.
"""
@with_kw mutable struct AgentParams <: AbstractAgentParams
    Z::Float64
    Idot_max_rel::Float64
    Idot_max_rel_emb::Float64
    X_emb_int::Float64
    H_p::Float64
    K_X::Float64
    
    """
    Initialize AgentParams from SpeciesParams `spc`.
    """
    function AgentParams(spc::AbstractParams)
        agn = new()
        agent_variability!(agn, spc)
        return agn
    end
end

"""
A `DEBParamCollection` contains global parameters `glb` and spc parameters `spc` (including TKTD-parameters). <br>
Initialize the default parameter collection with `DEBParamCollection()`.
"""
@with_kw mutable struct DEBParamCollection <: AbstractParamCollection
    glb::AbstractParams = GlobalParams()
    spc::AbstractParams = SpeciesParams()
    agn::Union{Nothing,AbstractParams} = nothing

    @assert length(spc.k_D_G) >= length(glb.C_W) "Length of k_D_G is not at least length of C_W in species"
    @assert length(spc.k_D_M) >= length(glb.C_W) "Length of k_D_M is not at least length of C_W in species"
    @assert length(spc.k_D_A) >= length(glb.C_W) "Length of k_D_A is not at least length of C_W in species"
    @assert length(spc.k_D_R) >= length(glb.C_W) "Length of k_D_R is not at least length of C_W in species"
    @assert length(spc.k_D_h) >= length(glb.C_W) "Length of k_D_h is not at least length of C_W in species"
    @assert length(spc.drc_functs_G) >= length(glb.C_W) "Length of drc_functs_G is not at least length of C_W in species"
    @assert length(spc.drc_functs_M) >= length(glb.C_W) "Length of drc_functs_M is not at least length of C_W in species"
    @assert length(spc.drc_functs_A) >= length(glb.C_W) "Length of drc_functs_A is not at least length of C_W in species"
    @assert length(spc.drc_functs_R) >= length(glb.C_W) "Length of drc_functs_R is not at least length of C_W in species"
    @assert length(spc.drc_functs_h) >= length(glb.C_W) "Length of drc_functs_h is not at least length of C_W in species"
    @assert length(spc.e_G) >= length(glb.C_W) "Length of e_G is not at least length of C_W in species"
    @assert length(spc.e_M) >= length(glb.C_W) "Length of e_M is not at least length of C_W in species"
    @assert length(spc.e_A) >= length(glb.C_W) "Length of e_A is not at least length of C_W in species"
    @assert length(spc.e_R) >= length(glb.C_W) "Length of e_R is not at least length of C_W in species"
    @assert length(spc.e_h) >= length(glb.C_W) "Length of e_h is not at least length of C_W in species"
    @assert length(spc.b_G) >= length(glb.C_W) "Length of b_G is not at least length of C_W in species"
    @assert length(spc.b_M) >= length(glb.C_W) "Length of b_M is not at least length of C_W in species"
    @assert length(spc.b_A) >= length(glb.C_W) "Length of b_A is not at least length of C_W in species"
    @assert length(spc.b_R) >= length(glb.C_W) "Length of b_R is not at least length of C_W in species"
    @assert length(spc.b_h) >= length(glb.C_W) "Length of b_h is not at least length of C_W in species"
end


"""
Definition of basic ABM object. <br>
Currently assumes that only a single species of type `AgentType` is simulated at a time.
"""
mutable struct DEBABM

    p::AbstractParamCollection # parameter collection
    t::Float64 # current simulation time
    dt::Float64 # timestep
    euler::Function # definition of the euler function for the given dt
    saveat::Float64
    
    u::ComponentVector # global state variables
    du::ComponentVector # global derivatives

    AgentID_count::Int64 # cumulative number of agents in the simulation
    agents::AbstractVector # population of agents
    aout::AbstractVector # agent output
    mout::AbstractVector # model output


    """
    Instantiate ABM from param collection `p`. 
    """
    function DEBABM(p::A; dt = 1/24, saveat = 1) where A <: AbstractParamCollection
        abm = new() # initialize ABM object
        abm.p = p # store the parameters
        abm.t = 0. # initialize simulation time
        abm.dt = dt # timestep is given as keyword argument
        abm.euler = defeuler(dt) # Euler function for the given dt
        abm.saveat = saveat
        
        abm.u = init_substates_global(p) # initialize the global substates
        abm.du = similar(abm.u) # initialize the global derivatives
        abm.du.X_p = 0.
        abm.du.C_W .= 0.
        
        abm.AgentID_count = 0 # set agent count
        initialize_agents!(abm) # initailize the population of agents
        abm.aout = [] # initialize agent output

        return abm
    end
end


@agent struct DEBAgent(GraphAgent)
    p::AbstractParamCollection
    u::ComponentVector
    du::ComponentVector
    AgentID::Int64
    dead::Bool
    cohort::Int64
    cum_repro::Int64

    """
        BaseAgent(p::AbstractParamCollection, AgentID::Int64)
    
    Initialize a base agent from parameters. 
    """
    function DEBAgent(abm::DEBABM)
        agent = new()

        agent.AgentID = abm.AgentID_count # assign AgentID
        abm.AgentID_count += 1 # increment AgentID counter

        agent.p = copy(abm.p) # parameters; TODO: #29 avoid copying everything
        agent.p.agn = AgentParams(agent.p.spc) # assign agent-level parameters (induces individual variability)
        initialize_statevars!(agent, abm) # initialize agent-level state variables
        agent.du = similar(agent.u) # initialize agent-level derivatives
        agent.du.glb = Ref(abm.du) # reference to global derivatives
        agent.du.agn = copy(agent.u.agn) # derivatives of agent substates have same shape as states

        agent.dead = false
        agent.cohort = 0
        agent.cum_repro += 1

        return agent
    end
end