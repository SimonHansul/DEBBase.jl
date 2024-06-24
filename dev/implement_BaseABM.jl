#=
# Base ABM implementation 

Implementation of the simplest possible generic version.<br>
Generic in terms of state variables and parameters, i.e. we don't need to make any chindes to this code if we add or remove variables or parameters
=#
@time "Loading packages" begin
    TAG = splitpath(@__FILE__)[end] |> x -> split(x, ".")[1] |> String    
    using Pkg; Pkg.activate("test")
    using Test

    using Plots, StatsPlots, Plots.Measures
    default(leg = false, lw = 1.5)
    
    using SHUtils
    using DataFrames, DataFramesMeta
    using ProgressMeter
    using Distributions
    using ComponentArrays
    using StaticArrays
    using OrdinaryDiffEq
    using Parameters, DEBParamStructs
    using Random
    using DEBFigures
    using Revise
    using DEBBase
    
    #include("../src/ModelFunctions.jl");
    #include("../src/Simulators.jl");
    #include("../src/IO.jl");
end

begin     
    p = DEBParamCollection()
    p.glb.N0 = 10

    p.glb.k_V = 0.5
    p.glb.Xdot_in = 3.
    p.glb.t_max = 30.
    p.spc.e_S = 0.5
    p.spc.b_S = 100.
    p.spc.eta_AR = 0.5

    @time mout, aout = ABM(p)

    pa = @df aout plot(
        groupedlineplot(
            :t, :S, :cohort, 
            ylabel = "S [μg C]", 
            leg = :topleft, label = hcat(unique(:cohort)...), 
            legendtitle = "Cohort", 
            legendtitlefontsize = 8)
    )

    pX = @df mout plot(:t, :X_p, ylabel = "Xₚ [μg C]")
    pN = @df mout plot(:t, :N_tot, ylabel = "N [#]")

    plot(pa, pX, pN, xlabel = "time [d]")
end   

# du.glb.X_p in dI!:  0x77cfa1ee1f628dc6
# du.glb.X_p in X_pdot_chemstat!: 0x9121e7fc286e9fbf


using Parameters, ComponentArrays

@with_kw mutable struct A
    x = ComponentVector(x1 = 1, x2 = 2)
end

mutable struct B
    x::ComponentVector
    
    function B(a::A)
        b = new()
        b.x = ComponentVector(
            xa = a.x,   
            xb = 1.
        )

        return b
    end
end


@with_kw mutable struct Model
    agents::AbstractVector = []
    X::Float64 = 1.
end

mutable struct Agent
    model::Model
end

m = Model()
a = Agent(m)


a.model.X += 1
m.X
