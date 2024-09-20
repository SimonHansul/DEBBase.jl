@with_kw mutable struct ABMParamCollection <: AbstractParamCollection
    glb::Union{AbstractParams,NamedTuple} = GlobalParams()
    spc::Union{AbstractParams,Vector{AbstractParams}} = [SpeciesParams()]
    agn::Union{Nothing,AbstractParams} = nothing
end

