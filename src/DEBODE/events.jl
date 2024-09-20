# events.jl
# functions used to construct callbacks callbacks in the ODE simulator 
# these can be used in the ABM be treating the agent as integrator

condition_juvenile(u, t, integrator) = u.X_emb # transition to juvenile when X_emb hits 0
function effect_juvenile!(integrator) 
    integrator.u.embryo = 0.
    integrator.u.juvenile = 1.
    integrator.u.adult = 0.
end
 
condition_adult(u, t, integrator) = u.H_p - u.H # condition to adult when H reaches H_p
function effect_adult!(integrator) 
    integrator.u.embryo = 0.
    integrator.u.juvenile = 0.
    integrator.u.adult = 1.
end
