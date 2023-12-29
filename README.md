# DEBBase.jl 

Functions to run reserveless DEB and DEB-TKTD models, including mixture effects. <br>
The idea of DEBBase is to serve as a simple base model for DEB-TKTD modelling. 
The purpose of DEBBase is explicitly not to provide functionality for parameter estimation or other kinds of analyses. This base model should be extendable and modular. 

## Features

- Custom datatypes to store DEB and global parameter sets
- Default parameter sets serve as a reference for implementation and testing
- Definition of ODE systems for DEB-TKTD with mixture exposure
- Functions to facilitate model input / output handling, e.g. converting ODE solution output to tidy dataframe

## Extending the model

You can use DEBBase to build your own extension of the base model. 
To do so, you first need to add `DEBBase.jl` and `Parameters` as a dependency to your project. <br>
In order to extend the DEBBase model, there are four implementation steps to take:


1. Define parameter objects
2. Define new model functions
3. Adapt the ODE system
4. Adapt the simulator function


In case you are adding new parameters to the model, 
you can define your own parameter struct like so: 
```Julia
using Parameters 
@with_kw mutable struct NewDEBParams <: AbstractParams
    Idot_max_rel_emb::Float64 = 22.9 # DEBBase parameter with default 
    k_E::Float64 = .1 # new parameter with default value variables in the model. this part ensures that mixture toxicity with arbitrary numbers of stressors can be simulated efficiently
    ...
end

@with_kw NewGlobalParams <: AbstractParams
    ...
end
```

`NewDEBParams` will replace `DEBBaseParams` and `NewGlobalParams` will replace `GlobalBaseParams`. So both need to include all the parameters, including those that were are already defined in `DEBBaseParams` and `GlobalBaseParams`. The definition of a new global parameter structure is optional, but a new DEB parameter structure should always be defined, even if no new parameters are introduced, because we will use multiple dispatch to differentiate between the base model and the extended model. <br>

Next, you need to define model functions which are additionally needed. For example, the function to calculate the derivative of some newly introduced state variable `E`:

```Julia
function Edot(
    deb::NewDEBParams;
    S::Float64
    )
    return deb.k_E * S
end
```

Some conventions for defining these functions in DEBBase:
- We assume dot notation for derivatives over time (`Edot`) 
- The only positional arguments are parameter objects and indices of stressors (matching the position of the stressor in the vector of exposure concentrations, vector of TKTD parameter values etc.)
- State variables which are arguments to the model function are given as keyword arguments with no default. <br>

Finally, you need to adapt the definition of the ODE system:

```
function NeWDEB!(du, u, p, t)
    glb, deb = p
    # unpack state variables
    X_p, X_emb, S, H, R, E = u

    S = max(0, S) # control for negative values

    life_stage = DEBBase.determine_life_stage(deb; H = H, X_emb = X_emb)
    
    idot = DEBBase.Idot(glb, deb; X_p = X_p, life_stage = life_stage, S = S)
    xdot = DEBBase.embryo(life_stage) ? 0. : glb.Xdot_in - idot
    xembdot = DEBBase.embryo(life_stage) ? -idot : 0.0
    adot = DEBBase.Adot(deb; Idot = idot)
    edot = Edot(deb; S = S)
    mdot = DEBBase.Mdot(deb; S = S)
    jdot = DEBBase.Jdot(deb; H = H)
    sdot = DEBBase.Sdot(deb; Adot = adot, Mdot = mdot) 
    hdot = DEBBase.Hdot(deb; Adot = adot, Jdot = jdot, adult = adult(life_stage))
    rdot = DEBBase.Rdot(deb; Adot = adot, Jdot = jdot, adult = adult(life_stage))

    # update du/dt
    du[01] = xdot
    du[02] = xembdot
    du[03] = sdot
    du[04] = hdot
    du[05] = rdot
    du[06] = edot
end
```

Compared to the original ODE system defined in `DEBBase.DEB!`, we made three modifications here: 

- Include `E` in the unpacking of state variables 
- Include `Edot` in the calculation of derivatives
- Include `edot` in the update of du/dt

In the definition of `NewDEB!`, we can easily see which variables are calculated according to the base model (`DEBBase.Idot`, `DEBBase.Adot`, etc.) and which functions are part of the extension (`Edot`). <br>


Finally, if new state variables have been introduced, the simulator also has to be adapted:

```Julia
function simulator(
    glb::AbstractParams,
    deb::NewDEBParams
    )
    
    u0 = [glb.Xdot_in, deb.X_emb_int, deb.X_emb_int * 0.01, 0., 0., 0.]
    tspan = (0, glb.t_max)
    prob = ODEProblem(NewDEB!, u0, tspan, (glb = glb, deb = deb))
    sol = solve(prob, reltol = 1e-6, abstol = 1e-10)
    simout = DataFrame(hcat(sol.t, hcat(sol.u...)'), [:t, :X, :X_emb, :S, :H, :R, :E])

    return simout
end
```

Here we have made the following changes:

- Added the initial value of the new state variable to `u0` and the corresponding symbol to the output data frame `simout`. 
- `ODEProblem` depends on `NewDEB!` instead of `DEB!`. <br>
- Exchanged `deb::AbstractParams` with `deb::NewDEBParams` in the function signature

The last step is optional - `AbstractParams` would also work, but specifiying the type is useful to add more methods to `simulator`.  
You might have many versions of `simulator` which all evaluate different parameter objects. 

### Tips on re-defining functions

A useful Julia function is `less()`, which prints the definition of a function. 
So you can for example use `less(DEBBase.DEB!)` to print the definition of the ODE system as given in the base model. 
Similarly, you can use `@less simulator(glb, deb)` to get the body for a specific method of `simulator` (depending on the type of `glb` and `deb`). <br>
Alternatively, most IDEs also have functions to navigate to the definition of a symbol (e.g. `Ctrl+Alt+Click` in VSCode - does not seem to always work well for functions defined in packages).

## TODO

- Add possibilities for simulating time-variable exposure (outsurce the detail work to a separate package?)