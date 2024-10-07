skipinf(x::Array) = x[isfinite.(x)]

vectify(x) = parse.(Float64, split(split(split(x, "[")[end], "]")[1]," "))

function which_in(x::AbstractString, possibilities::Vector{String}; none_found_return_val="") 
    idxs = findall(z->occursin(z, x), possibilities)
    if length(idxs)>0
        return possibilities[idxs[1]]
    else
        return none_found_return_val
    end
end


"""
    geomrange(v::Vector{Float64}; length=50)

Geometric series created from range of values within a vector.
"""
geomrange(v::Vector{Float64}; length=50) = 10 .^ range(log10(minimum(v)), log10(maximum(v)), length=length)

"""
    geomrange(a::Real, b::Real; length=50)

Geometric series created from two extreme values.
"""
geomrange(a::Real, b::Real; length=50) = 10 .^ range(log10(a), log10(b); length=length)

"""
    diffvec(x::AbstractVector)

Calculate difference along a vector, inserting NaN as first element (analogous to Pandas' diff method).
"""
diffvec(x::AbstractVector) = vcat([NaN], diff(x))

"""
    fround(x; sigdigits=2)

Formatted rounding to significant digits (omitting decimal point when appropriate). 
Returns rounded number as string.
"""
function fround(x; sigdigits=2)
    xround = string(round(x, sigdigits=sigdigits))
    if xround[end-1:end]==".0"
        xround = string(xround[1:end-2])
    end
    return xround
end


"""
    drop_missing(df::AbstractDataFrame; verbose=false)::DataFrame

Drop all rows with missing values from data frame.
"""
function drop_missing(df::AbstractDataFrame; verbose=false)::DataFrame
    n0 = nrow(df)
    df2 = df[completecases(df),:]
    dn = nrow(df2)-n0
    if verbose
        @info("Dropped $dn of $n0 rows containing missing values.")
    end
    return df2
end

function replace_missing!(
    df::AbstractDataFrame,
    cols::Vector{Symbol};
    replace_val = 0.0
    )
    for col in cols
        df[ismissing.(df[:,col]),col] .= replace_val
    end
    return df
end

function get_treatment_names(exposure::Vector{Vector{N}}, stressor_names::Array{Symbol,1}) where N <: Number
    treatment_type = ["co"]
    treatment_level = [0]
    treatment = ["co"]

    treatment_level_counter = 0
    stressor_pre = "co"
    for (i,x) in enumerate(exposure[2:end])
        sum(x.>0) == 1 ? stressor = string(stressor_names[x.>0][1]) : stressor = "mix"
        if stressor == stressor_pre
            treatment_level_counter += 1
        else
            treatment_level_counter = 1
        end
        push!(treatment_type, string(stressor))
        push!(treatment_level, treatment_level_counter)
        push!(treatment, stressor * string(treatment_level_counter))
        stressor_pre = stressor
    end
    return treatment_type, treatment_level, treatment
end

function lab(v::Vector{R}; kwargs...) where R <: Real
    return hcat(unique(fround.(v; kwargs...))...)
end

function wrappend(file::String, data::DataFrame, step::Int64)
    if (isfile(file)==false)&(step>1)
        error("Attempt to append to non-existing file: step>1 but file does not exist.")
    end
    if step == 1
        CSV.write(file, data, append=false)
    else
        CSV.write(file, data, append=true)
    end
end

function read_W3C(file_path::AbstractString; kwargs...)::DataFrame
    meta = []
    core_data = []

    # Open the file for reading
    open(file_path, "r") do file
        for line in eachline(file)
            # Check if the line starts with a hashtag
            if startswith(line, "#")
                # Process metadata
                push!(meta, line)
            else
                # Process core data
                push!(core_data, line)
            end
        end
    end

    core_data_str = join(core_data, "\n")
    core_data_table = CSV.File(IOBuffer(core_data_str); kwargs...) |> DataFrame

    return core_data_table
end

function ismin(x::AbstractVector)
    return x .== minimum(x)
end

get_pkg_version(pkg, deps) = OrderedDict(deps).vals |> x ->filter(d -> d.name == pkg, x) |> x -> x[1].version
