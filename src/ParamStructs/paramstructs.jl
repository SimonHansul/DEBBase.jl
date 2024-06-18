#=
# Abstract parameter types
=#

abstract type AbstractParams end
abstract type AbstractParamCollection end

abstract type AbstractSpeciesParams <: AbstractParams end
abstract type AbstractGlobalParams <: AbstractParams end

copy(theta::AbstractParams) = typeof(theta)([getproperty(theta, x) for x in fieldnames(typeof(theta))]...)
copy(theta::AbstractParamCollection) = typeof(theta)([getproperty(theta, x) for x in fieldnames(typeof(theta))]...)