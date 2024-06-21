#=
# Functions to generate parameter structs from base parameter structs

These functions make it possible so that we can take an existing parameter struct (e.g. the DataType `DEBBaseParams`) 
and define a new mutable struct with additional parameters. 
The only function which is relevant for users here is `defparams()`.
=#


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
    # FIXME: Figure out how "<: abstracttype" has to be included in the quote...
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
Example:

    defparams(:Xenopus, DEBBaseParams, (:gamma => 0.5,)) # define a parameter struct Xenopus, with additional parameter :gamma.
    deb = Xenopus() # create an instance of Xenopus with default parameters

The additional parameter can in turn have an arbitrary type.
"""
function defparams(Name, BaseParams, newparams)
    eval(genstruct(Name, BaseParams, newparams))
end
