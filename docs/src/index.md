# DEBBase.jl: A basis for Dynamic Energy Budget modelling in Julia


DEBBase includes a number of submodules. 

The core of the package is formed by `DEBBase.DEBODE` and `DEBBase.ABC`.

`DEBODE` defines a baseline DEB-TKTD model and provides the core functionality to run simulations. 
`ABC` contains code for parameter estimation using Sequential Monte Carlo Approximate Bayesian Computation. 


## DEBODE: DEB-TKTD model as system of ordinary differential equations 

```@autodocs
Modules = [DEBBase.DEBODE]
```

## ABC: Parameter inference using Sequential Monte Carlo Approximate Bayesian Computation
