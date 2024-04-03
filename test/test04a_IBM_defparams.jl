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
    using Parameters, SpeciesParamstructs
end

@test begin
    plt = plot(
        layout = (1,2), 
        leftmargin = 2.5mm, bottommargin = 2.5mm, 
        xlabelfontsize = 10, 
        size = (600,350),
        xlabel = "Time since fertilization (d)",
        ylabel = ["Structure" "Reproduction buffer"]
        )

    yhat_fix = DEBBase.simulator(DEBParamCollection())
    @df yhat_fix plot!(plt, :t, :S, subplot = 1, color = :black, lw = 2)
    @df yhat_fix plot!(plt, :t, :R, subplot = 2, color = :black, lw = 2)
    
    hyperZ = Truncated(Normal(1., 0.1), 0, Inf)
    yhat_var = @replicates DEBBase.simulator(DEBParamCollection(spc = SpeciesParams(Z = hyperZ))) 10
    @df yhat_var plot!(plt, :t, :S, group = :replicate, alpha = .25, subplot = 1, c = :viridis)
    @df yhat_var plot!(plt, :t, :R, group = :replicate, alpha = .25, subplot = 2, c = :viridis) 

    true
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
abstract type AbstractAgent end

"""
DEBBase Agent. <br>
Each agent owns a reference to its associated parameter collection.
"""
mutable struct BaseAgent <: AbstractAgent
    pcmn::Base.RefValue{DEBParamCollection} # reference to the common parameter collection
    pown::ComponentVector
    u::ComponentVector
    du::ComponentVector
    active::Bool

    function BaseAgent(pcmn::Ref{A}) where A <: AbstractParamCollection # generate agent from reference to paramter collection
        a = new()
        a.pcmn = pcmn
        a.pown = DEBBase.initialize_pown()
        DEBBase.agent_variability!(a.pown, a.pcmn)
        a.u = DEBBase.initialize_statevars(a.pcmn, a.pown)
        a.du = similar(a.u)
        a.active = false
        return a
    end
end

function activate!(a::BaseAgent)
    a.active = true
end

function deactivate!(a::BaseAgent)
    a.active = false
end

@test begin
    Θ = DEBParamCollection()
    Θ_ref = Ref(Θ)
    a = BaseAgent(Θ_ref)    

    true
end




#=
## Implementing a model object 


=#

abstract type AbstractABM end


θ = DEBParamCollection()
θ_ref = Ref(θ)

agents = Vector{AbstractAgent}(undef, 10)
agents[1] = BaseAgent(θ_ref)

"""
Definition of basic ABM object. <br>
Currently assumes that only a single species of type `AgentType` is simulated at a time.
"""
@with_kw mutable struct ABM <: AbstractABM
    theta::AbstractParamCollection # Parameter collection
    agents # Agents
    t::Float64

    """
    Instantiate ABM from param collection `theta`. 
    """
    function ABM(theta::A; AgentType = BaseAgent, N_max = Int(1e4)) where A <: AbstractParamCollection
        abm = new() # initialize ABM object
        abm.theta = theta # 
        abm.t = 0.

        abm.agents = Vector{AbstractAgent}(undef, N_max)
        
        for i in eachindex(abm.agents)
            abm.agents[i] = AgentType(Ref(theta))
            if i <= theta.glb.N0
                activate!(abm.agents[i])
            end
        end

        return abm
    end
end


abm = ABM(θ)
println(abm)