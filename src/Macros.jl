
if occursin("terminal", abspath(PROGRAM_FILE))
    @info("Running setup")
    using Pkg; Pkg.activate("test")
    using DEBParamStructs
    using Parameters
    using DataFrames, DataFramesMeta
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
Macro to generate
"""
macro defparams(Name, field, val)
    eval(genstruct(Name, field, val))
end


@defparams PASF gamma 0.5
@defparams wholenotherparams blabla [x -> 1 + x, x -> 3x]

PASF().gamma


wholenotherparams()