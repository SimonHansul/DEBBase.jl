# DEBBase.jl: A basis for Dynamic Energy Budget modelling in Julia


DEBBase is a basic package to perform Dynamic Energy Budget Toxicokinetic-Toxicodynamic (DEB-TKTD) modelling in Julia. 
It includes a number of submodules:

The core of the package is formed by `DEBBase.DEBODE` and `DEBBase.ABC`.

`DEBODE` defines a baseline DEB-TKTD model and provides the core functionality to run simulations. 
`ABC` contains code for parameter estimation using Sequential Monte Carlo Approximate Bayesian Computation. 


## DEBODE: DEB-TKTD model as system of ordinary differential equations 

```@autodocs
Modules = [DEBBase.DEBODE]
```

## ABC: Parameter inference using Sequential Monte Carlo Approximate Bayesian Computation


```@autodocs
Modules = [DEBBase.ABC]
``` 


## DoseResponse: A collection of Dose-Response functions

```@autodocs
Modules = [DEBBase.DoseResponse]
``` 


## Figures: Recipes and utilities for plotting model output

```@autodocs
Modules = [DEBBase.Figures]
``` 

## ParamStructs: Infrastructure to organize parameter structures

By default, DEBBase uses mutable structs to store parameters. <br>


```@autodocs
Modules = [DEBBase.ParamStructs]
``` 

## Utils: Other utility functions

```@autodocs
Modules = [DEBBase.Utils]
``` 





