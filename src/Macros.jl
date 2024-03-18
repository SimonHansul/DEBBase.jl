
if occursin("terminal", abspath(PROGRAM_FILE))
    @info("Running setup")
    using Pkg; Pkg.activate("test")
    using DEBParamStructs
    using Parameters
    using DataFrames, DataFramesMeta
    @time using DEBBase
end

include("Structures.jl")

"""
    genstruct(Name, field, val)

This function generates an expression which in turn generates a parameter structure. \n
This makes up the core of the `@defparams` macro, which is simply a wrapper around `genstruct()`.
"""
function genstruct(Name, field, val)
    quote 
        @with_kw mutable struct $Name
            $field = $val
        end
    end
end


"""
Macro to generate parameter structure
"""
macro defparams(Name, field, val)
    eval(genstruct(Name, field, val))
end

@defparams SomeNewParams gamma 0.5
SomeNewParams().gamma

#=
Ok, we know now how to define a mutable struct inlcuding defaults within a macro. 
The next step is to extend this to multiple fields.
=#

function genstruct(Name, fields, vals)
    exprs = Expr[]
    for (field, val) in zip(fields, vals)
        push!(exprs, :( $(Symbol(field)) = $val))
    end

    quote 
        @with_kw mutable struct $Name
            $(exprs...)
        end
    end
end

"""
Macro to generate parameter structure
"""
macro defparams(Name, fields, vals)
    eval(genstruct(Name, fields, vals))
end

expr = genstruct(:Bloahr, [:k, :g], [1., 3.])

obj = eval(expr)

obj()


(exprs...)


@defparams MultipleParams [:kappa, :gamma] [0.5, 0.8]


@with_kw mutable struct Ehm
    x = 1; y = 2
end


Ehm()

fields = [:p1, :p2]
vals = [0.5, 1.0]

e = Expr(:tuple, [:( $field = $val ) for (field, val) in zip(fields, vals)]...)

eval(e)





wholenotherparams()