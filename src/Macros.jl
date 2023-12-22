"""
Unpack vector of state variables. 
Global Parameters have to contain `statevars` (a `Vector{Symbol}` indicating the names of state variables) 
and `sttypes` (a `Vector{DataType}` indicating the types of state variables).
$(TYPEDSIGNATURES)
"""
macro unpack(u, glb::AbstractParams)
    quote
        for (i, u_i) in enumerate(u)
            eval(:($(glb.statevars[i]) = $(glb.sttypes[i])($(u_i))))
        end
    end
end