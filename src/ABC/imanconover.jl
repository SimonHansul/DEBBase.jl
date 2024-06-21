
"""
Induce rank correlations to independent samples, using Iman-Conover algorithm. 
This function was adopted from the MCHammer package. 
---
**MCHammer was released with the following License statement: **\n
Copyright 2019-2021, Technology Partnerz Ltd. and Eric Torkia \n
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: \n
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. \n
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. \n
"""
function iman_conover(ar, cor_mat)
    n_trials = size(ar)[1]
    if typeof(ar) == Array{Float64,2}
          ar = DataFrame(ar)
    end

    #Define how many columns of ISNs are required
    array_dims = size(ar,2)

    #calc cholesky transform to transfer correlation to ISNs
    P = cholesky(cor_mat)

    #Normal() Returns Standard Normals (ISNs)
    R = rand(Normal(),n_trials,array_dims)
    ISN_Matrix = R*P.U
    ISN_Matrix_DF = DataFrame(ISN_Matrix, :auto)

    #apply ranks to create independant correlation rankings matrix
    ISN_Ranked = []
    for i = 1:array_dims
          temp_ranks = ordinalrank(ISN_Matrix_DF[!,i])
          push!(ISN_Ranked, temp_ranks)
    end

    #Reindex the array of samples using the ISN_Ranks. Sort(Array)[OrderingVector]
    final_array=[]
    for i = 1:array_dims
          sorted_array = sort(ar[:,i])[ISN_Ranked[i]]
          push!(final_array, sorted_array)
    end
    return hcat(final_array...)
end
