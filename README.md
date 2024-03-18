# DEBBase.jl 

Functions to run reserveless DEB and DEB-TKTD models, including mixture effects. <br>
The purpose of DEBBase is to serve as a simple base model for DEB-TKTD modelling and provide 
the basic infrastucture to perform DEB-TKTD modelling. 
The purpose of DEBBase is explicitly not to provide functionality for parameter estimation or other kinds of analyses. 

## Features

- Custom data types to store DEB and global parameter sets
- Default parameter sets serve as a reference for implementation and testing
- Definition of ODE systems for DEB-TKTD with mixture exposure
    - Reserveless DEB 
    - TK with account for body size
    - Flexible TD: DRC functions are parameters, not hard-coded (two-parameter log-logistic by default)
    - Arbitrary number of chemical stressors
    - Mixture toxicity simulated based on Independent Action

- Functions to facilitate model input / output handling, e.g. converting ODE solution output to tidy data frame

## Quickstart

### Installation

DEBBase can be installed via Github, but the unregistered dependencies have to be installed manually first:

```Julia#
# load the package manager
using Pkg

# activate project in current directory 
# (or whatever project you want to add DEBBase to)
Pkg.activate(".") 

# add unregistered dependencies
Pkg.add(url = "https://github.com/SimonHansul/DEBParamStructs.jl")
Pkg.add(url = "https://github.com/SimonHansul/DoseResponse.jl")

# add DEBBase
Pkg.add(url = "https://github.com/SimonHansul/DEBBase.jl")
```

### Usage

The following Code will simulate the DEB model based on the given default parameters:

```Julia
using DEBBase
out = DEBBase.simulator(BaseParams())
```

A parameter collection (`BaseParamCollection`) contains two sets of parameters: 
- `glb::GlobalBaseParams`: These are global parameters, such as the simulated timespan, food input rate, etc.
- `deb::DEBBaseParams`: The DEB and TKTD parameters.

---

# DEB Modelling: languages & tools (I)

How does DEBBase.jl fit into the landscape of DEB-TKTD/IBM tools? Here is an attempt to briefly describe the strenghts of different tools:


| Tool            | Language | Open source | Individual-level | Population-level | TKTD | Mixtures |Published | Performance |Parameter estimation |
|-----------------|----------|-------------|------------------|------------------|------|----------|----------|-------------|--------|
| Debtool         | Matlab   | -^***       | +                | -                |  -   |    -/?   |+         | ?           |+       |
| DEBBase         | Julia    | +           | +                | +                |  +   |     +    |-         | +           |+* | 
| MEMpy           | Python   | +           | +                | +                |  +   |    +     |-         | -           |-** |
| Netlogo DEB-IBM | Netlogo  | +           | +/-              | +                | +    |    +     |+         | -           |-** |   

* Added through DEBABC package

** These tools can of course be used for parameter inference, but there are no convenient built-in functions

*** DEBtool itself is open source, but Matlab is not


Or, to summarize:

| Tool            | Unique features                                      |
|-----------------|------------------------------------------------------|
| Debtool         | Tried and tested manifestation of standard DEB model |
| DEBBase         | Extensible performant implementation with mixture TKTD & integrated individual-based modelling capabilities |
| MEMpy           | Easy-to-use Python implementation for individuals and populations |
| Netlogo DEB-IBM | Probably the most used DEB-IBM implementation |




## TODO

- Add possibilities for simulating time-variable exposure (outsurce the detail work to a separate package?)
