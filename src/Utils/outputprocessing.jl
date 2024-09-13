

"""
    sol_to_mat(sol::O)::Matrix{Float64}
Convert ODE solution object to matrix.
"""
function sol_to_mat(sol::O)::Matrix{Float64} where O <: ODESolution
    return hcat(sol.t, hcat(sol.u...)')
end

"""
Convert ODE solution object to output data frame.
"""
function sol_to_df(sol::O)::DataFrame where O <: ODESolution
    simout = DataFrame(
        sol_to_mat(sol), # ODE output converted to matrix
        extract_colnames(sol) # column names inferred from component array keys
    )
    return simout
end

"""
Calculate relative response, given scalar values.
"""
function rr(x::R, x_ref::R)::Float64 where R <: Real
    if x_ref == 0.
        return 1.
    else
        return x / x_ref
    end
end

function rr(x::R, x_ref::Missing)::Missing where R <: Real
    return missing
end

function robustmean(x)
    xfilt = filter(xi -> isfinite(xi), x)

    if length(xfilt)==0
        return NaN
    end

    return mean(xfilt)
end

"""
Calculate the relative responses. \n
Positional arguments: 
- `sim::AbstractDataFrame`: results
- `response_vars::Vector{Symbol}`: response variables for which to calculate the relative responses
- `treatment_var::Symbol`: Column indicating the treatment. Column values can be numerical or categorical, but `identify_control` kwarg has to be specified in the latter case

Keyword arguments: 
- `groupby_vars::Vector{Symbol}`: relative response will be conditioned on these variables (e.g. time, separate experiments...). Empty by default.
- `identify_control`: function used to identify reference values from `treatment_var`. By default, this is `minimum()` (assuming that column values in `treatment_var` are numerical).
"""
function relative_response(
    sim::D, 
    response_vars::Vector{Symbol},
    treatment_var::Symbol; 
    groupby_vars::Vector{Symbol} = Symbol[],
    identify_control = minimum
    ) where D <: AbstractDataFrame

    #=
    Calculation of the conditional control mean
    =#
    reference = sim[sim[:,treatment_var] .== identify_control(sim[:,treatment_var]),:] |> # extract the control
    x -> groupby(x, groupby_vars) |> # group by conditioning variables
    x -> combine(x) do df 
        refvals = DataFrame() # (sub-)dataframe of reference values
        for var in response_vars # iterate over response values
            var_ref = Symbol(String(var) * "_ref") # get the reference column name
            refvals[!,var_ref] = [robustmean(df[:,var])] # calculate the conditional control mean
        end
        return refvals
    end

    #=
    Calculation of the relative response
    =#
    sim = leftjoin(sim, reference, on = groupby_vars) # add references as new columns
    for var in response_vars # for every response variable
        y_var = Symbol("y_" * String(var)) # get the relative response column name
        var_ref = Symbol(String(var) * "_ref") # get the reference column name
        sim[!,y_var] = [rr(row[var], row[var_ref]) for row in eachrow(sim)] # calculate the relative response
        select!(sim, Not(var_ref)) # drop the reference column
    end
    return sim
end

"""
    idcol!(sim::AbstractDataFrame, col::Symbol, val)

Add identifier column with name `col` and value `val` to results DataFrame `sim`
"""
function idcol!(sim::AbstractDataFrame, col::Symbol, val)
    sim[!,col] .= val
end

"""
    idcol!(sim::Any, col::Symbol, val)

Add identifier column to a results object `sim`, assuming that `sim` is some collection of `DataFrame`s (e.g. Tuple, Vector).
"""
function idcol!(sim::Any, col::Symbol, val)
    for df in sim
        idcol!(df, col, val)
    end
end
