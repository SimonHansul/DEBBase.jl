# DEBBase.jl: A basis for Dynamic Energy Budget modelling in Julia

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://simonhansul.github.io/DEBBase.jl)

This is a Julia package to define and run DEB-TKTD models, including mixture effects. <br>
The purpose of DEBBase is to provide a simple base model for DEB-TKTD modelling and provide 
the basic infrastucture to perform DEB-TKTD modelling. 
The purpose of DEBBase is explicitly not to provide functionality for parameter estimation or other kinds of analyses. 

## Key features

- Custom data types to store DEB and global parameter sets
- Default parameter sets serve as a reference for implementation and testing
- Definition of ODE systems for DEB-TKTD with mixture exposure
    - Reserveless DEB 
    - Minimal toxikokinetics
    - Flexible toxicodynamics: Dose-response functions can be swapped out via the species-specific parameters.
    - Mixture toxicity simulated based on Independent Action or Damage Addition
    - Simulation of an arbitrary number of mixture exposures with arbitrary PMoA combinations
    - Functions to facilitate model input / output handling, e.g. converting ODE solution output to tidy data frame
    - Likelihood-free Bayesian parameter inference using Sequential Monte Carlo Approximate Bayesian Computation

## Installation & Quickstart

DEBBase is currently not registered, but can be installed via Github. 

```Julia
using Pkg; Pkg.add("https:/github.com/simonhansul/debbase.jl")
```

The following Code will simulate the DEB model based on the given default parameters:

```Julia
using DEBBase
params = DEBParamCollection()
out = DEBBase.simulator(DEBParamCollection())
```

`params` will be a nested mutable struct containing a set of default parameters. <br>
The only purpose of the default parameters is to provide a place to get started and a reference during development.<br>
Check the [documentation](docs/build/index.html) for details on all implemented functions. 
(Note that within the documentation of each function, you can click on `source` to inspect the corresponding source code.) 
You can also inspect the test scripts for more simulation examples.

## Earmarked future features (road to v1.0.0)

- Simulation of time-variable exposure
- Using callbacks for life stage transitions
- Integration of the individual-based implementation of the base model
- Integration with ModellingToolkit.jl

## Acknowledgements 

This package was developed as part of the research project AmphiDEB, funded by the European Food Safety Authority.