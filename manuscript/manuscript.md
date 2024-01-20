# The DEBBase.jl ecosystem: accessible and extensible simulation of life-history, population dynamics and stressor effects using Dynamic Energy Budgets

Simon Hansul, Andreas Focks

## Context

Models based on Dynamic Energy Budget theory are gaining increasing attention in environmental sciences, especially in combination with Toxicokientic-Toxicodynamic (TKTD) models. DEB-TKTD models allow to extract information about the mechanisms of toxicity from life-history data, and to subsequently extrapolate the effects of chemical stressors from effects on the indivdual level to effects on population and community dynamics. As such, they are both valuable for basic scientific research and are of increasing interest for the application in environmental risk assessment. <br>

## Motivation

Several implementations of DEB and DEB-TKTD models exist. 
However, there is need for a accessible yet flexible implementation 
of DEB-TKTD models, which allows to efficiently simulate individual life-history 
as well as population dynamics. Furthermore, models should be modifiable within a clear framework,  which should also enhance proper documentation of the modifications. The DEBBase.jl ecosystem aims to provide such an implementation. <br> 

## Design principles

DEBBase.jl is designed to be lean and modular. 
That means, instead of providing as much functionality as possible in a single package, 
DEBBase.jl only provides the definition of a base model, default parameters and simulators for the base model. Simulations are either be executed by passing the model on to an ODE solver, 
or by passing them on to an individual-based simulation framework. This is a novel feature of DEBBase.jl compared to existing DEB-TKTD model implementations, since those focus either on simulating individual life-history (e.g. DEBtool, BYOM) or on simulating population dynamics (e.g. the Netlogo implementation by Martin et al. (2013)). <br>

Parameter estimation is outsourced to the separate DEBABC.jl, but other Julia packages could be used to perform parameter estimation (e.g. Optim.jl). Vice versa, DEBABC.jl can also be used to perform parameter estimation for models which are not part of DEBBase.jl.
DEBBase.jl can be used as a dependency to implement more specific DEB models, as we will later demonstrate using a DEB model for amphibians.

## Simulating life-history

With DEBBase.jl installed in the Julia environment, we can conduct a first simulation with just two lines of code:
```Julia
using DEBBase
simout = simulator(BaseParamCollection())
```

Here, the first line loads the DEBBase package. The second line calls the `simulator` function, which takes an object of type `BaseParamCollection` as argument. 
This will assume a set of default parameters if no parameters are specified, as done above. <br>
The definition of `BaseParamCollections` in turn contains two fields:

- `glb`: An instance of type `GlobalBaseParams` which defines global parameters. These include quantities such as the timespan to simulate, nutrient input rates, and external chemical concentrations. The global parameters in combination define the environmental or experimental scenario to be simulated. The default values simulate life-history of the default DEB organism at *ad libitum* feeding conditions and without exposure to chemical stressors over a timespan of 21 days.
- `deb`: An instance of type `DEBBaseParams`, which defines the DEB and TKTD parameters. The default organism roughly approximates the growth and reproduction of *Daphnia magna*.

Figure 1 shows the growth and reproduction returned by the default parameter collection.

![img](../plots/basetest_growthrepo.png) <br>
*Figure 1: Structure $S$ and reproduction $R$ of the default organism in the default scenario.*

The default parameters can be modified by specifying exactly the parameters which should be modified, e.g:

```Julia
p = BaseParamCollection()
p.deb.kappa = 0.8
simout = simulator(p)
```

This simulates the default parameters, except that the DEB parameter $\kappa$ (indicating the fraction of assimilated nutrients allocated to growth and somatic maintenance) is changed to 0.8 Alternatively, the following code achieves the same:

```Julia
simout = simulator(BaseParamCollection(deb = DEBBaseParams(kappa = 0.8)))
```

The same applies to global parameters, e.g. to set the food input rate:
```Julia
simout = simulator(BaseParamCollection(glb = GlobalBaesParams(Xdot_in = 800.)))
```

## Simulating effects of chemical stressors

The TKTD component is implemented so that mixtures with an arbitrary number of mixture components can be simulated. For this reason, TKTD parameters are stored as Vectors, where each element of the Vector represents a mixture component. Furthermore, we assume that each mixture component can act via arbitrary combinations of Physiological Modes of Action (PMoA), represented by a separate parameter Vector for each PMoA. The PMoAs are indicated by the suffixes G = decrease in growth efficiency, M = increase in maintenance costs, A = decrease in assimilation efficiency, R = decrease in reproduction efficiency, h = lethal effects according to GUTS-RED-SD. <br>
In the case of the TK parameters, each element of the Vector ` is a scalar value, 
because we assume a simple TK model with a single parameter, the dominant rate konstant $k_{D,z,j}$, specific for chemical stressor $z$ and PMoA $j$.

```

```

## Parameter estimation using DEBABC.jl


## Simulating population dynamics


## Extending the base model


## Accessibility


## Outlook

- Environmental factors, e.g. temperature
- Variable exposure concentrations
- Mixtures: Damage addition, combinations of DA and IA, interaction factors
- Sensitivity analysis

