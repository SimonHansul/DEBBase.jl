"""
Create shared y-labels for plot with grid layout, placing labels only on the left-most subplots.
"""
function gridylabel(label::String, numrows::Int64, numcols::Int64)::Matrix{String}
    labelstring = [label] # the label as vector
    numemptystrings = numcols-1 # the number of empty strings needed per row
    emptystrings = repeat([""], numemptystrings) # vector of empty strings
    labels = repeat(vcat(labelstring, emptystrings), numrows) # repeat for each row
    labels = hcat(labels...) # convert to matrix
    return labels
end

function gridxlabel(label::String, numrows::Int64, numcols::Int64)::Matrix{String}
    labelstring = [label] # the label as vector
    numemptyrows = numrows-1 # the number of rows with empty strings
    emptystrings = repeat([""], numemptyrows * numcols) # vector of empty strings
    lastrow = repeat(labelstring, numcols) # the last row, containing the actual labels
    labels = vcat(emptystrings, lastrow) # repeat for each row
    labels = hcat(labels...) # convert to matrix
    return labels
end

"""
Save plot plt with name, assuming known global TAG and plotsdir.
"""
function saveplt(plt, name)
    savefig(plt, plotsdir("$(TAG)_$(name).pdf"))
    savefig(plt, plotsdir("$(TAG)_$(name).png"))
end