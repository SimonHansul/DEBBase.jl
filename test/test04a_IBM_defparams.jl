@time "Loading packages" begin
    using Pkg; Pkg.activate("test")
    using Test
    using Plots, StatsPlots, Plots.Measures
    using SHUtils
    using DataFrames, DataFramesMeta
    using ProgressMeter
    using Distributions
    default(leg = false, lw = 1.5)
    TAG = splitpath(@__FILE__)[end] |> x -> split(x, ".")[1] |> String    

    using Revise
    using DEBBase
    using Parameters, DEBParamStructs
end

begin
    plt = plot(
        layout = (1,2), 
        leftmargin = 2.5mm, bottommargin = 2.5mm, 
        xlabelfontsize = 10, 
        size = (600,350),
        xlabel = "Time since fertilization (d)",
        ylabel = ["Structure" "Reproduction buffer"]
        )

    y = DEBBase.simulator(BaseParamCollection())
    @df y plot!(plt, :t, :S, subplot = 1, color = :black, lw = 2)
    @df y plot!(plt, :t, :R, subplot = 2, color = :black, lw = 2)
    
    hyperZ = Truncated(Normal(1., 0.1), 0, Inf)
    y = @replicates DEBBase.simulator(BaseParamCollection(deb = DEBBaseParams(Z = hyperZ))) 10
    @df y plot!(plt, :t, :S, group = :replicate, alpha = .25, subplot = 1, c = :viridis)
    @df y plot!(plt, :t, :R, group = :replicate, alpha = .25, subplot = 2, c = :viridis) 

    display(plt)
end


#=
## Implementing an Agent object. 

What does an agent need?

- A reference to common paramters `pcmn`
- Agen-specific parameters `pown`
- State variables `u`
- Derivatives `du`

- A function to initialize the agent
=#

using ComponentArrays


"""
DEBBase Agent. <br>
Each agent owns a reference to its associated parameter collection.
"""
mutable struct BaseAgent
    pcmn::Base.RefValue{BaseParamCollection} # reference to the common parameter collection
    pown::ComponentVector
    u::ComponentVector
    du::ComponentVector

    function BaseAgent(pcmn::Ref{A}) where A <: AbstractParamCollection # generate agent from reference to paramter collection
        a = new()
        a.pcmn = pcmn
        a.pown = DEBBase.initialize_pown()
        DEBBase.agent_variability!(a.pown, a.pcmn)
        a.u = DEBBase.initialize_statevars(a.pcmn, a.pown)
        a.du = similar(a.u)
        return a
    end
end


theta = BaseParamCollection()
thtref = Ref(theta)

a = BaseAgent(thtref)

DEBBase.DEB!()


macro expand(innerfunc, outerfunc)
    quote
        $outerfunc(du, u, p, t) = begin
            $(esc(innerfunc))(du, u, p, t)
        end
    end
end

function X_pindot!(du, u, p, t)
    du.X_p = p[1].x.glb.Xdot_in
end

expmodel = @expand DEBBase.DEB! X_pindot!

@macroexpand expmodel

dump(DEBBase.DEB!)

