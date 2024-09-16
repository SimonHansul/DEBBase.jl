using Pkg; Pkg.activate("docs")
using Documenter
using DEBBase

makedocs(
    sitename = "DEBBase.jl", 
    modules = [DEBBase, DEBBase.ABC, DEBBase.DEBODE],
    format = Documenter.HTML()
    )