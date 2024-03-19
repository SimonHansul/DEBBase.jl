
if occursin("terminal", abspath(PROGRAM_FILE))
    @info("Running setup")
    using Pkg; Pkg.activate("test")
    using DEBParamStructs
    using Parameters
    using DataFrames, DataFramesMeta
    @time using DEBBase
end

# TODO: this should probably live in DEBParamStructs.jl?

include("Structures.jl")


"""
Generate parameter fields with defaults as Vector of expressions.
"""
function paramexpr(params)
    exprs = Expr[]
    for pair in params
        push!(exprs, :( $(Symbol(pair.first)) = $(pair.second)))
    end
    return exprs
end

function paramexpr(params::DataType)
    instance = params()
    exprs = Expr[]

    for param in fieldnames(params)
        defval = getproperty(instance, param)
        expr = :( $(Symbol(param)) = $(defval))
        push!(exprs, expr)
    end

    return exprs
end


"""
    genstruct(Name, BaseParams::DataType, newparams::Tuple{Pair{Symbol,Any},...})

Generate an expression to define a parameter structure named `Name`, given base paramater struct `BaseParams` and added parameters `newparams` (the latter is a tuple of pairs).
Note that `defparams` is a wrapper around `genstruct` which evaluates the expression.
"""
function genstruct(Name::Symbol, BaseParams::DataType, newparams; abstracttype = AbstractParams)
    
    exprs = vcat(
        paramexpr(BaseParams),
        paramexpr(newparams)
    )

    quote 
        @with_kw mutable struct $Name
            $(exprs...)
        end
    end
end


"""
    defparams(Name::Symbol, BaseParams::DataType, newparams::Tuple{Pair{Symbol,Any},...})
Define a parameter structure named `Name`, given base paramater struct `BaseParams` and added parameters `newparams` (the latter is a tuple of pairs).
"""
function defparams(Name, BaseParams, newparams)
    eval(genstruct(Name, BaseParams, newparams))
end



defparams(:Lymnaea, DEBBaseParams, (:M_L => 0.1,))
