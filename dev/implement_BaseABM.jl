#=
# Base ABM implementation 

Implementation of the simplest possible generic version.<br>
Generic in terms of state variables and parameters, i.e. we don't need to make any chagnes to this code if we add or remove variables or parameters
=#


@time "Loading packages" begin
    using Pkg; Pkg.activate("test")
    TAG = splitpath(@__FILE__)[end] |> x -> split(x, ".")[1] |> String    
    using Test

    using Plots, StatsPlots, Plots.Measures
    default(leg = false, lw = 1.5)
    
    using SHUtils
    using DataFrames, DataFramesMeta
    using ProgressMeter
    using Distributions
    using ComponentArrays
    using StaticArrays
    using Parameters, DEBParamStructs
    using Random
    using DEBFigures
    using Revise

    include(("../src/Params.jl"))
    #include("../src/ModelFunctions.jl");
    #include("../src/Simulators.jl");
    #include("../src/IO.jl");
end


"""
Macro to initialize parameter structure with Id = 1. 
Usage: 
```
@agent struct Params(NoSpaceAgent)
    x = 1
end

p = @init Params
```
"""
macro init(type)
    quote
        $type(id = 1)
    end
end

