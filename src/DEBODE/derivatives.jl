"""
Clip negative values.
"""
function clipneg(x::Float64)
    return sig(x, 0., 0., x)
end

"""
Sigmoid switch function. 
`y_left` and `y_right` are the function values left and right of the threshold `x_thr`.

"""
@inline function sig(
    x::Float64, 
    x_thr::Float64,
    y_left::Float64, 
    y_right::Float64; 
    beta::Float64 = 1e16
    )::Float64

    return 1 / (1 + exp(-beta*(x - x_thr))) * (y_right - y_left) + y_left
end

@inline function functional_response(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Float64

    let X_V = u.X_p / p.glb.V_patch # convert food abundance to concentration
        return X_V / (X_V + p.ind.K_X) # calculate type II functional response
    end
end

"""
Calculate ingestion rate. 
Embryos (X_emb <= 0) take up resources from the vitellus X_emb. 
Juveniles and adults (X_emb > 0) feed on the external resource X_pcmn.

"""
@inline function Idot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Nothing

    u.f_X = functional_response(du, u, p, t)
    
    du.I_emb = sig( # uptake from vitellus
        u.X_emb, # uptake from vitellus depends on mass of vitellus
        0., # the switch occurs when vitellus is used up 
        0., # when the vitellus is used up, there is no uptake
        (Complex(u.S)^(2/3)).re * p.ind.Idot_max_rel; # when the vitellus is not used up, uptake from vitellus occurs
        beta = 1e20 # for switches around 0, we need very high beta values
        )

    du.I_p = sig( # uptake from external resource p
        u.X_emb, # ingestion from external resource depends on mass of vitellus
        0., # the switch occurs when the vitellus is used up  
        u.f_X * p.ind.Idot_max_rel * (Complex(u.S)^(2/3)).re, # when the vitellus is used up, ingestion from the external resource occurs
        0.; # while there is still vitellus left, there is no uptake from the external resource
        beta = 1e20 # again we have a switch around 0, requiring very high beta
        )

    du.I = du.I_emb + du.I_p

    du.X_p -= du.I_p
    du.X_emb = -du.I_emb

    return nothing
end

"""
Assimilation flux
"""
@inline function Adot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real)::Nothing

    du.A = du.I * p.spc.eta_IA * u.y_A

    return nothing
end

"""
Somatic maintenance flux
"""
@inline function Mdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real)::Nothing

    du.M = u.S * p.spc.k_M * u.y_M

    return nothing
end

"""
Maturity maintenance flux

"""
@inline function Jdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Nothing

    du.J = u.H * p.spc.k_J * u.y_M

    return nothing
end

"""
Positive somatic growth
"""
@inline function Sdot_positive(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Float64

    return p.spc.eta_AS * u.y_G * (p.spc.kappa * du.A - du.M)
end

"""
Negative somatic growth
"""
@inline function Sdot_negative(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Float64 

    return -(du.M / p.spc.eta_SA - p.spc.kappa * du.A)
end

function Sdot(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Float64 

    return sig(
        p.spc.kappa * du.A, # growth depends on maintenance coverage
        du.M, # switch occurs based on maintenance costs
        Sdot_negative(du, u, p, t), # left of the threshold == maintenance costs cannot be covered == negative growth
        Sdot_positive(du, u, p, t) # right of the threshold == maintenance costs can be covered == positive growth
    )
end

"""
Somatic growth, including application of shrinking equation.
"""
@inline function Sdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Nothing

    du.S = Sdot(du, u, p, t)

    return nothing
end

"""
    S_max_histdot!(
        du::ComponentArray,
        u::ComponentArray,
        p::AbstractParamCollection,
        t::Real
        )::Nothing

Reference structure `S_max_hist`, which is used as a measure of energetic state.

This is the amount of structure an individual of the given age has under ideal conditions, 
i.e. `f = 1` and `y_z = 1`, with the exception of `y_G`. 
Values of `y_G != 1` are included in the calculation of `S_max_hist`, so that a slowing down of growth 
"""
@inline function S_max_histdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Nothing

    #du.S_max_hist = p.eta_AS * u.y_G * (DEBBase.sig(u.X_emb, 0., p.ind.Idot_max_rel, p.ind.Idot_max_rel_emb; beta = 1e20) * Complex(u.S ^(2/3)).re - p.spc.k_M * u.S_max_hist)

    u.S_max_hist = sig(
        u.S,
        u.S_max_hist,
        u.S_max_hist,
        u.S;
    )

    return nothing
end

"""
Maturation flux. 

Maturity is dissipated energy and can therefore not be burned to cover maintenance costs. 

Currently there are no consequences for an organism not covering maturity maintenance. 
If this turns out to be an issue, we might consider to add a damage pool ``D_H``,
where the amount of maturity maintenance that could not be covered is accumulated. 
This might then lead to a fitness penalty depending on ``D_H``, for example in the form of additional 
mortality or embryonic hazard. 

Such rules are however likely species-specific and should be evaluated in the light of a more precise problem definition.
"""
@inline function Hdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Nothing 

    du.H = sig(
        u.H, # maturation depends on maturity
        p.ind.H_p, # switch occurs at maturity at puberty H_p
        clipneg(((1 - p.spc.kappa) * du.A) - du.J), # maturation for embryos and juveniles
        0., # maturation for adults
    )

    return nothing
end

"""
Update the current estimate of H_b. 
The current estimate of maturity at birth is equal to current maturity for embryos, 
and will be fixed to the current value upon completion of embryonic development.
This way, we obtain H_b as an implied trait and can use it later (for example in `abj()`).
"""
@inline function H_bdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Nothing

    du.H_b = sig(
        u.X_emb, # estimate depends on embryonic state
        0., # switch occurs when vitellus is gone
        0., # post-embryonic: H_b stays fixed
        du.H # embryonic: H_b tracks H
    )

    return nothing
end

"""
Reproduction flux.
"""
@inline function Rdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Nothing 

    du.R = sig(
        u.H, # reproduction depends on maturity
        p.ind.H_p, # switch occurs at maturity at puberty H_p
        0., # reproduction for embryos and juveniles
        clipneg(u.y_R * p.spc.eta_AR * ((1 - p.spc.kappa) * du.A - du.J)) # reproduction for adults
    )

    return nothing
end

"""
    Qdot!(
        du::ComponentArray,
        u::ComponentArray,
        p::AbstractParamCollection,
        t::Real
        )::Nothing 

Calculation of the total dissipation flux, equal to the sum of maintenance costs and overheads paid for assimilation, mobilization, maturation, growth and reproduction.
"""
function Qdot!(
    du::ComponentArray,
    u::ComponentArray,
    p::AbstractParamCollection,
    t::Real
    )::Nothing 

    # dissipation fluxes for the individual processes
    let Qdot_A, Qdot_S, Qdot_R
        
        Qdot_A = du.I * (1 - p.spc.eta_IA)
        Qdot_S = du.S >= 0 ? du.S * (1 - p.spc.eta_AS) / p.spc.eta_AS : du.S * (p.spc.eta_SA - 1)
        Qdot_R = du.R * (1 - p.spc.eta_AR) / p.spc.eta_AR
        
        du.Q =  Qdot_A + Qdot_S + Qdot_R + du.M + du.J + du.H
    end

    return nothing
end

"""
    C_Wdot!(
        du::ComponentVector,
        u::ComponentVector,
        p::AbstractParamCollection,
        t::Real
        )::Nothing 

Change in external concentrations. 
Currently simply returns zeros because time-variable exposure is not yet implemented. 
# TODO: extend to account for time-variable exposure.
"""
@inline function C_Wdot!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )::Nothing 

    du.C_W = zeros(length(u.C_W)) # constant exposure : derivative is 0

    return nothing
end

"""
Change in environmental resource abundance, simulating a chemostat.
"""
function X_pdot!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )::Nothing 

    du.X_p = p.glb.Xdot_in - p.glb.k_V * u.X_p

    return nothing
end

"""
Constant concentration of external chemical stressor:
"""
function C_Wdot_const!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )::Nothing

    du.C_W = zeros(length(u.C_W))

    return nothing
end


"""
TK for spc-TKTD model, including effect of surface area to volume ratio and dilution by growth. 
If `D` and is given as a Vector, TK is only stressor-specific but not PMoA-specific. 
"""
@inline function Ddot!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )::Nothing 

    for z in eachindex(u.C_W)
        @inbounds du.D_G[z] = sig(u.X_emb, 0., p.spc.k_D_G[z] * (u.C_W[z] - u.D_G[z]), 0.)
        @inbounds du.D_M[z] = sig(u.X_emb, 0., p.spc.k_D_M[z] * (u.C_W[z] - u.D_M[z]), 0.)
        @inbounds du.D_A[z] = sig(u.X_emb, 0., p.spc.k_D_A[z] * (u.C_W[z] - u.D_A[z]), 0.)
        @inbounds du.D_R[z] = sig(u.X_emb, 0., p.spc.k_D_R[z] * (u.C_W[z] - u.D_R[z]), 0.)
        @inbounds du.D_h[z] = sig(u.X_emb, 0., p.spc.k_D_h[z] * (u.C_W[z] - u.D_h[z]), 0.)
    end

    return nothing
end

"""
Response to chemical stressors, assuming independent action for mixtures.
"""
@inline function y_z!(
    du::ComponentVector,
    u::ComponentVector,
    p::AbstractParamCollection,
    t::Real
    )::Nothing 

    @inbounds u.y_G = prod([p.spc.drc_functs_G[z](u.D_G[z], (p.spc.e_G[z], p.spc.b_G[z])) for z in 1:length(u.C_W)]) # combined relative responses for sublethal effects per PMoA
    @inbounds u.y_M = prod([p.spc.drc_functs_M[z](u.D_M[z], (p.spc.e_M[z], p.spc.b_M[z])) for z in 1:length(u.C_W)])
    @inbounds u.y_A = prod([p.spc.drc_functs_A[z](u.D_A[z], (p.spc.e_A[z], p.spc.b_A[z])) for z in 1:length(u.C_W)])
    @inbounds u.y_R = prod([p.spc.drc_functs_R[z](u.D_R[z], (p.spc.e_R[z], p.spc.b_R[z])) for z in 1:length(u.C_W)])

    @inbounds u.h_z = sum([p.spc.drc_functs_h[z](u.D_h[z], (p.spc.e_h[z], p.spc.b_h[z])) for z in 1:length(u.C_W)]) # hazard rate

    return nothing
end



#"""
#    Hbj(H::Float64, X_emb::Float64, H_b::Float64, H_j::Float64, p_b::Float64, p_j::Float64)::Float64
#
#Mautrity-driven metabolic acceleration from birth to maturity threshold `H_j` (metamorphosis). 
#This form of metabolic acceleration assumes that metabolic rate constants change as maturation progresses 
#(e.g. due to changes in gene expression). 
#
#We assume that some baseline parameter `p` has value `p_b` at birth and `p_j` at metamorphosis.
#Between birth and metamorphosis, the current value of `p` is the maturity-weighted mean of `p_b` and `p_j`.
#"""
#function Hbj(H::Float64, X_emb::Float64, H_b::Float64, H_j::Float64, p_b::Float64, p_j::Float64)::Float64
#    w_b = (H_j - H) / (H_j - H_b) # weight for p_b
#    w_j = 1 - w_b # weight for p_j
#    p_bj = mean([p_b, p_j], Weights([w_b, w_j])) # p_bj, i.e. value between birth and maturity
#    
#    p = DEBBase.sig( # post-metamorphosis: value stays constant at p_j
#        H,
#        H_j,
#        DEBBase.sig( # embryonic: value stays constant at p_b
#            X_emb, 
#            0., 
#            p_bj, 
#            p_b),
#        p_j
#    )
#
#    return p
#end

<<<<<<< HEAD:src/ModelFunctions.jl
    return nothing
end

"""
Mautrity-driven metabolic acceleration from birth to maturity threshold `H_j` (metamorphosis). 
We assume that some baseline parameter `p` has value `p_b` at birth and `p_j` at metamorphosis.
Between birth and metamorphosis, the current value of `p` is the maturity-weighted mean of `p_b` and `p_j`.
"""
function Hbj(H::Float64, X_emb::Float64, H_b::Float64, H_j::Float64, p_b::Float64, p_j::Float64)::Float64
    w_b = (H_j - H) / (H_j - H_b) # weight for p_b
    w_j = 1 - w_b # weight for p_j
    p_bj = mean([p_b, p_j], Weights([w_b, w_j])) # p_bj, i.e. value between birth and maturity
    
    p = DEBBase.sig( # post-metamorphosis: value stays constant at p_j
        H,
        H_j,
        DEBBase.sig( # embryonic: value stays constant at p_b
            X_emb, 
            0., 
            p_bj, 
            p_b),
        p_j
    )

    return p
end


"""
    reproduce_opportunistic!(agent::AbstractAgent, abm::AbstractABM)::Nothing

Reproduction according to an opportunistic reproduction strategy. 
Offspring is created whenever there is enough energy in the reproduction buffer.
"""
function reproduce_opportunistic!(agent::AbstractAgent, abm::AbstractABM)::Nothing
    
    let num_offspring = trunc(agent.u.agn.R / agent.p.agn.X_emb_int) # calculate the number of offspring produced
        for _ in 1:num_offspring # for the given number of offspring agents
            offspring = abm.p.glb.AgentType(abm) # initialize a new agent
            offspring.u.agn.cohort = agent.u.agn.cohort + 1 # the offspring agent is part of a new cohort
            push!(abm.agents, offspring) # add offspring to agents vector
            agent.u.agn.R -= agent.p.agn.X_emb_int # remove energy from reproduction buffer
        end
        agent.u.agn.cum_repro += num_offspring
    end

    return nothing
end

"""
    ingestion!(agent::AbstractAgent, abm::AbstractABMagent)::Nothing

Update resource abundance in model `abm` based on ingestion flux of `agent`.
"""
function ingestion!(agent::AbstractAgent, abm::AbstractABM)::Nothing
    
    abm.u.X_p = max(0, abm.u.X_p - agent.du.agn.I_p * abm.dt)
    
    return nothing
end

"""
    stochastic_death!(agent::AbstractAgent, abm::AbstractABM, h_x::Float64)::Nothing

Apply stochastic death model to agent, with respect to hazard rate `h_x`. 

Records the cause of death encoded in a number, given by keyword argument `causeofdeath`.
"""
function stochastic_death!(agent::AbstractAgent, abm::AbstractABM, h_x::Float64; causeofdeath = 0.)::Nothing
    if rand() >= exp(-h_x / abm.dt)
        agent.u.agn.dead = 1.
        agent.u.agn.causeofdeath = causeofdeath
    end

    return nothing
end

function die!(agent::AbstractAgent, abm::AbstractABM)::Nothing
    stochastic_death!(agent, abm, agent.u.agn.h_S; causeofdeath = 1.) # starvation mortality from loss of structure
    return nothing
end

function N_tot!(abm::AbstractABM)::Nothing
    abm.u.N_tot = length(abm.agents)
    return nothing
end
=======
>>>>>>> 68ffa3666dd1b56a58a769b08e648a1011cd23be:src/DEBODE/derivatives.jl

"""
    DEBODE!(du, u, p, t)

Definition of base model as a system of ordinary differential equations. 
This model definition is suitable for simulating the life-history of a single organism in conjecture with DifferentialEquations.jl.
"""
function DEBODE!(du, u, p, t)::Nothing

    #### physiological responses
    y_z!(du, u, p, t) # response to chemical stressors
<<<<<<< HEAD:src/ModelFunctions.jl
    h_S!(du, u, p, t) # starvation mortality
=======
>>>>>>> 68ffa3666dd1b56a58a769b08e648a1011cd23be:src/DEBODE/derivatives.jl

    #### auxiliary state variables (record cumulative values)
    Idot!(du, u, p, t)
    Adot!(du, u, p, t) 
    Mdot!(du, u, p, t) 
    Jdot!(du, u, p, t)
    Qdot!(du, u, p, t)

    #### major state variables
    Sdot!(du, u, p, t) # structures
    S_max_histdot!(du, u, p, t) # reference structure
    Hdot!(du, u, p, t) # maturity 
    H_bdot!(du, u, p, t) # estimate of maturity at birth
    Rdot!(du, u, p, t) # reproduction buffer
    X_pdot!(du, u, p, t) # resource abundance
    Ddot!(du, u, p, t) # damage
    C_Wdot!(du, u, p, t) # external stressor concentration 

    return nothing
end