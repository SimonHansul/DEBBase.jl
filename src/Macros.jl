"""
Unpack vector of state variables
"""
macro unpack(u, glb::AbstractParams)
    quote
        for (i, u_i) in enumerate(u)
            eval(:($(glb.statevars[i]) = $(glb.sttypes[i])($(u_i))))
        end
    end
end