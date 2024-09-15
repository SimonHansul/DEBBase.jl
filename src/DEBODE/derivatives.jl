#derivatives.jl
# A collection of derivative functions which is used to compose pre-defined models, and can be re-used for custom models.
# Note that the ODE system itself is treated as a species parameter `odefunc`, where the derivatives are given as a Vector of Functions, in the order in which they are called.

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
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Float64

    let X_V = u.X_p / p.glb.V_patch # convert food abundance to concentration
        return X_V / (X_V + u.K_X) # calculate type II functional response
    end
end

"""
    dI!(
        du::ComponentArray,
        u::ComponentArray,
        p::Union{AbstractParamCollection,NamedTuple},
        t::Real
        )::Nothing
        
Calculate ingestion rate. 
Embryos (X_emb <= 0) take up resources from the vitellus X_emb. 
Juveniles and adults (X_emb > 0) feed on the external resource X_pcmn.
"""
@inline function dI!(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Nothing

    u.f_X = functional_response(du, u, p, t)
    
    du.I_emb = u.embryo * (Complex(u.S)^(2/3)).re * u.Idot_max_rel_emb # uptake from vitellus for embryos
    du.I_p = (1 - u.embryo) * u.f_X * u.Idot_max_rel * (Complex(u.S)^(2/3)).re # uptake from external resource for all other life stages
    du.I = du.I_emb + du.I_p # total uptake is the sum (though normally it is either one or the other)

    du.X_p -= du.I_p # change in external resource abundance
    du.X_emb = -du.I_emb # change in vitellus

    return nothing
end

"""
Assimilation flux
"""
@inline function dA!(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real)::Nothing

    du.A = du.I * u.eta_IA

    return nothing
end

"""
Somatic maintenance flux
"""
@inline function dM!(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real)::Nothing

    du.M = u.S * u.k_M

    return nothing
end

"""
Maturity maintenance flux

"""
@inline function dJ!(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Nothing

    du.J = u.H * u.k_J

    return nothing
end

"""
Positive somatic growth.
"""
@inline function Sdot_positive(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Float64

    return u.eta_AS * (u.kappa * du.A - du.M)
end

"""
Negative somatic growth, assuming that the residual ``\\kappa``-assimiliation flux and structure will be combined to pay somatic maintenance.
"""
@inline function Sdot_negative(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Float64 

    return -(du.M / p.spc.eta_SA - u.kappa * du.A)
end


"""
Somatic growth, accounting for the possibility of shrinking.
"""
function Sdot(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Float64 

    return sig(
        u.kappa * du.A, # growth depends on maintenance coverage
        du.M, # switch occurs based on maintenance costs
        Sdot_negative(du, u, p, t), # left of the threshold == maintenance costs cannot be covered == negative growth
        Sdot_positive(du, u, p, t) # right of the threshold == maintenance costs can be covered == positive growth
    )
end

"""
Somatic growth, including application of shrinking equation.
"""
@inline function dS!(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Nothing

    du.S = Sdot(du, u, p, t)

    return nothing
end

"""
    dS_max_hist!(
        du::ComponentArray,
        u::ComponentArray,
        p::Union{AbstractParamCollection,NamedTuple},
        t::Real
        )::Nothing

Reference structure `S_max_hist`, which is used as a measure of energetic state.

This is the amount of structure an individual of the given age has under ideal conditions, 
i.e. `f = 1` and `y_z = 1`, with the exception of `y_G`. 
Values of `y_G != 1` are included in the calculation of `S_max_hist`, so that a slowing down of growth 
"""
@inline function dS_max_hist!(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Nothing

    
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
@inline function dH!(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Nothing 

    du.H = (1 - u.adult) * clipneg(((1 - u.kappa) * du.A) - du.J)

    return nothing
end

"""
Update the current estimate of H_b. 
The current estimate of maturity at birth is equal to current maturity for embryos, 
and will be fixed to the current value upon completion of embryonic development.
This way, we can obtain H_b as an implied trait for use during the simulation.
"""
@inline function dH_b!(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
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
@inline function dR!(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Nothing 

    du.R = u.adult * clipneg(u.y_R * u.eta_AR * ((1 - u.kappa) * du.A - du.J)) # reproduction for adults

    return nothing
end

"""
    dQ!(
        du::ComponentArray,
        u::ComponentArray,
        p::Union{AbstractParamCollection,NamedTuple},
        t::Real
        )::Nothing 

Calculation of the total dissipation flux, equal to the sum of maintenance costs and overheads paid for assimilation, mobilization, maturation, growth and reproduction.
"""
function dQ!(
    du::ComponentArray,
    u::ComponentArray,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Nothing 

    # dissipation fluxes for the individual processes
    let Qdot_A, Qdot_S, Qdot_R
        
        Qdot_A = du.I * (1 - u.eta_IA)
        Qdot_S = du.S >= 0 ? du.S * (1 - u.eta_AS) / u.eta_AS : du.S * (p.spc.eta_SA - 1)
        Qdot_R = du.R * (1 - u.eta_AR) / u.eta_AR
        
        du.Q =  Qdot_A + Qdot_S + Qdot_R + du.M + du.J + du.H
    end

    return nothing
end

"""
    dC_W!(
        du::ComponentVector,
        u::ComponentVector,
        p::Union{AbstractParamCollection,NamedTuple},
        t::Real
        )::Nothing 

Change in external concentrations. 
Currently simply returns zeros because time-variable exposure is not yet implemented. 
# TODO: extend to account for time-variable exposure.
"""
@inline function dC_W!(
    du::ComponentVector,
    u::ComponentVector,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Nothing 

    du.C_W = zeros(length(u.C_W)) # constant exposure : derivative is 0

    return nothing
end

"""
Change in environmental resource abundance, simulating a chemostat.
"""
function dX_p!(
    du::ComponentVector,
    u::ComponentVector,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Nothing 

    du.X_p = p.glb.Xdot_in - p.glb.k_V * u.X_p

    return nothing
end


"""
TK for spc-TKTD model, including effect of surface area to volume ratio and dilution by growth. 
If `D` and is given as a Vector, TK is only stressor-specific but not PMoA-specific. 
"""
@inline function dD!(
    du::ComponentVector,
    u::ComponentVector,
    p::Union{AbstractParamCollection,NamedTuple},
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
@inline function y_z_IndependentAction!(
    du::ComponentVector,
    u::ComponentVector,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Nothing 

    @inbounds u.y_G = prod([p.spc.drc_functs_G[z](u.D_G[z], (p.spc.e_G[z], p.spc.b_G[z])) for z in 1:length(u.C_W)]) # combined relative responses for sublethal effects per PMoA
    @inbounds u.y_M = prod([p.spc.drc_functs_M[z](u.D_M[z], (p.spc.e_M[z], p.spc.b_M[z])) for z in 1:length(u.C_W)])
    @inbounds u.y_A = prod([p.spc.drc_functs_A[z](u.D_A[z], (p.spc.e_A[z], p.spc.b_A[z])) for z in 1:length(u.C_W)])
    @inbounds u.y_R = prod([p.spc.drc_functs_R[z](u.D_R[z], (p.spc.e_R[z], p.spc.b_R[z])) for z in 1:length(u.C_W)])
    @inbounds u.h_z = sum([p.spc.drc_functs_h[z](u.D_h[z], (p.spc.e_h[z], p.spc.b_h[z])) for z in 1:length(u.C_W)]) # hazard rate

    return nothing
end


"""
Response to chemical stressors, assuming damage addition for mixtures. <br>

The PMoA-specific damage values for each stressor are added up,

``D_j = \\sum_{z=1}^{n}{D_{j,z}}``

, and the response is caluclated assuming identical DRC parameters for all stressors. 
Currently, there are no additional weight factors implemented (assuming that the model is fitted to single-substance data only).

"""
function y_Z_DamageAddition!(
    du::ComponentVector,
    u::ComponentVector,
    p::Union{AbstractParamCollection,NamedTuple},
    t::Real
    )::Nothing 
    
    u.y_G = p.spc.drc_functs_G(sum(u.D_G), (p.spc.e_G[z], p.spc.b_G[1]))
    u.y_M = p.spc.drc_functs_G(sum(u.D_G), (p.spc.e_G[z], p.spc.b_G[1]))
    u.y_A = p.spc.drc_functs_G(sum(u.D_G), (p.spc.e_G[z], p.spc.b_G[1]))
    u.y_R = p.spc.drc_functs_G(sum(u.D_G), (p.spc.e_G[z], p.spc.b_G[1]))
    u.h_z = p.spc.drc_functs_G(sum(u.D_G), (p.spc.e_G[z], p.spc.b_G[1]))

    return nothing
end


"""
    Hbj(H::Float64, X_emb::Float64, H_b::Float64, H_j::Float64, p_b::Float64, p_j::Float64)::Float64

Mautrity-driven metabolic acceleration from birth to maturity threshold `H_j` (metamorphosis). 
We assume that some baseline parameter `p` has value `p_b` at birth and `p_j` at metamorphosis.
Between birth and metamorphosis, the current value of `p` is the maturity-weighted mean of `p_b` and `p_j`. 
Not used in the base model.
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
Calculate Arrhenius temperature correction factor.
"""
function tempcorr!(
    du::ComponentVector, 
    u::ComponentVector, 
    p::Union{AbstractParamCollection,NamedTuple}, 
    t::Real)::Nothing
    
    u.y_T = exp(p.spc.T_ref/u.T * ((p.spc.T_ref - p.spc.T_A)/(u.T - p.spc.T_A)))
    
    return nothing
end


"""
Apply stressors to baseline parameter values. 
Temperature correction affects rate parameters, 
other stressors affect parameters according to their respective PMoA.
"""
function apply_stressors!(du, u, p, t)

    u.eta_AS = p.spc.eta_AS_0 * u.y_G 
    u.k_M = p.spc.k_M_0 * u.y_M * u.y_T
    u.k_J = p.spc.k_J_0 * u.y_M * u.y_T
    u.eta_IA = p.spc.eta_IA_0 * u.y_A
    u.eta_AR = p.spc.eta_AR_0 * u.y_R

    u.Idot_max_rel_emb = p.agn.Idot_max_rel_emb_0 * u.y_T
end


"""
    DEBBase!(du, u, p, t) # putting the model together

The DEBBase model derivatives. <br>

This model is used the default by DEBBase, and can be regarded as a version of DEBkiss with maturity (cf. Tjalling Jager's DEBkiss book).

Feeding is simulated explicitly and all state variables, including structure, have the same dimension. 
This means that the density of dry mass does not appear in the model - structure is a mass or energy pool like everyting else.
The implication is that the surface area-specific ingestion rate as the dimension ``dim(u)^{1/3} d^{-1}``, 
where ``dim(u)`` is the chosen model currency (e.g. dry mass or energy).

The model implementation allows for simulations of mixtures with an arbitrary number of substances, 
assuming independent action. The default dose-response is a log-logistic function, which is ``1-ln``-tansformed for maintenance costs 
(increasing function with lower limit at 1).
"""
function DEBBase!(du, u, p, t) # putting the model together

    y_z_IndependentAction!(du, u, p, t) # calculate response to chemical stressors
    tempcorr!(du, u, p, t) # calculate response to 
    apply_stressors!(du, u, p, t) # apply stressors to baseline parameters

    #### auxiliary state variables (record cumulative values)
    dI!(du, u, p, t)
    dA!(du, u, p, t) 
    dM!(du, u, p, t) 
    dJ!(du, u, p, t)
    dQ!(du, u, p, t)

    #### major state variables
    dS!(du, u, p, t)
    dS_max_hist!(du, u, p, t)
    dH!(du, u, p, t)
    dH_b!(du, u, p, t)
    dR!(du, u, p, t)
    dX_p!(du, u, p, t)
    dD!(du, u, p, t)
    dC_W!(du, u, p, t)
end