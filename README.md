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
out = simulator(BaseParamCollection())
```

A parameter collection (`BaseParamCollection`) contains two sets of parameters: 
- `glb::GlobalBaseParams`: These are global parameters, such as the simulated timespan, food input rate, etc.
- `deb::DEBBaseParams`: The DEB and TKTD parameters.

The command `BaseParamCollection()` initiates a parameter collection with the default parameters, which approximate the life history of *D. magna* (model currency $\mu g\ C$). Correspondingly, `GlobalBarams()` and `DEBBaseParams()` initiates only the default `glb`  and `deb` parameters, respectively. <br>
Parmeters can be modified two ways. Either by passing them on as keyword arguments:
```Julia
# initialize default DEB parameters, but change kappa
deb = DEBBaseParams(kappa = 0.5)
theta = BaseParamCollection(deb = deb) # initialize parameter collection with altered DEB parameters
```

Or, by accessing the parameter entry of a previously initialized parameter structure:
```Julia
deb = DEBBaseParams() # initialize default DEB parameters
deb.kappa = 0.5 # change kappa
setproperty!(deb, :kappa, 0.5) # this does the same as the line above
theta = BaseParamCollection(deb = deb) # initialize parameter collection
```


## TODO

- Add possibilities for simulating time-variable exposure (outsurce the detail work to a separate package?)
