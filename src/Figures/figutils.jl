"""
Create shared y-labels for plot with grid layout, placing labels only on the left-most subplots.
"""
function gridylabel(label::AbstractString, numrows::Int64, numcols::Int64)::Matrix{String}
    labelstring = [label] # the label as vector
    numemptystrings = numcols-1 # the number of empty strings needed per row
    emptystrings = repeat([""], numemptystrings) # vector of empty strings
    labels = repeat(vcat(labelstring, emptystrings), numrows) # repeat for each row
    labels = hcat(labels...) # convert to matrix
    return labels
end

function gridxlabel(label::AbstractString, numrows::Int64, numcols::Int64)::Matrix{String}
    labelstring = [label] # the label as vector
    numemptyrows = numrows-1 # the number of rows with empty strings
    emptystrings = repeat([""], numemptyrows * numcols) # vector of empty strings
    lastrow = repeat(labelstring, numcols) # the last row, containing the actual labels
    labels = vcat(emptystrings, lastrow) # repeat for each row
    labels = hcat(labels...) # convert to matrix
    return labels
end

"""
Assign some attribute value `attr` to the left-most column of a multi-plot and 
another attribute value `nullattr` to  all other subplots.

Example:

```Julia
using Plots, Plots.Measures
# sets the left left column margin to 2.5mm and the other ones to 0mm
plot(plot(), plot(), plot(), leftmargin = gridyattr(2.5mm, 0mm, 1, 3), layout = (1,3)) 
```

"""
function gridyattr(attr, nullattr, numrows::Int64, numcols::Int64)
    attrvec = [attr] # the label as vector
    numnullattrs = numcols-1 # the number of empty strings needed per row
    nullattrs = repeat([nullattr], numnullattrs) # vector of empty strings
    attrs = repeat(vcat(attrvec, nullattrs), numrows) # repeat for each row
    attrs = hcat(attrs...) # convert to matrix
    return attrs
end

"""
    gridxattr(attr, nullattr, numrows::Int64, numcols::Int64)::Matrix{String}

Analogous to `gridyattr`.
"""
function gridxattr(attr, nullattr, numrows::Int64, numcols::Int64)
    attrvec = [attr] # the label as vector
    numnullattrs = numrows-1 # the number of rows with empty strings
    nullattrs = repeat([nullattr], numnullattrs * numcols) # vector of empty strings
    lastrow = repeat(attrvec, numcols) # the last row, containing the actual labels
    attrs = vcat(nullattrs, lastrow) # repeat for each row
    attrs = hcat(attrs...) # convert to matrix
    return attrs
end


"""
Save plot plt with name, assuming known global TAG and plotsdir.
"""
function saveplt(plt, name)
    savefig(plt, plotsdir("$(TAG)_$(name).pdf"))
    savefig(plt, plotsdir("$(TAG)_$(name).png"))
end