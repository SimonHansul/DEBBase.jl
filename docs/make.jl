using Pkg; Pkg.activate("docs")
using Documenter
using DEBBase
using DEBBase.ABC

makedocs(
    sitename = "DEBBase: Dynamic Energy Budget modelling in Julia", 
    modules = [DEBBase, DEBBase.ABC, DEBBase.DEBODE],
    format = Documenter.HTML()
    )