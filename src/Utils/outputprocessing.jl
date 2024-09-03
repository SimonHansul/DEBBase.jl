

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


"""
Calculate the relative responses. \n
Positional arguments: 
- `fit::AbstractDataFrame`: results
- `response_vars::Vector{Symbol}`: response variables for which to calculate the relative responses
- `treatment_var::Symbol`: Column indicating the treatment. Column values can be numerical or categorical, but `identify_control` kwarg has to be specified in the latter case

Keyword arguments: 
- `groupby_vars::Vector{Symbol}`: relative response will be conditioned on these variables (e.g. time, separate experiments...). Empty by default.
- `identify_control`: function used to identify reference values from `treatment_var`. By default, this is `minimum()` (assuming that column values in `treatment_var` are numerical).
"""
function relative_response(
    fit::D, 
    response_vars::Vector{Symbol},
    treatment_var::Symbol; 
    groupby_vars::Vector{Symbol} = Symbol[],
    identify_control = minimum
    ) where D <: AbstractDataFrame

    #=
    Calculation of the conditional control mean
    =#
    reference = fit[fit[:,treatment_var] .== identify_control(fit[:,treatment_var]),:] |> # extract the control
    x -> groupby(x, groupby_vars) |> # group by conditioning variables
    x -> combine(x) do df 
        refvals = DataFrame() # (sub-)dataframe of reference values
        for var in response_vars # iterate over response values
            var_ref = Symbol(String(var) * "_ref") # get the reference column name
            refvals[!,var_ref] = [mean(df[:,var])] # calculate the conditional control mean
        end
        return refvals
    end

    #=
    Calculation of the relative response
    =#
    fit = leftjoin(fit, reference, on = groupby_vars) # add references as new columns
    for var in response_vars # for every response variable
        y_var = Symbol("y_" * String(var)) # get the relative response column name
        var_ref = Symbol(String(var) * "_ref") # get the reference column name
        fit[!,y_var] = [rr(row[var], row[var_ref]) for row in eachrow(fit)] # calculate the relative response
        select!(fit, Not(var_ref)) # drop the reference column
    end
    return fit
end

"""
    idcol!(Yhat::AbstractDataFrame, col::Symbol, val)

Add identifier column with name `col` and value `val` to results DataFrame `Yhat`
"""
function idcol!(Yhat::AbstractDataFrame, col::Symbol, val)
    Yhat[!,col] .= val
end

"""
    idcol!(Yhat::Any, col::Symbol, val)

Add identifier column to a results object `Yhat`, assuming that `Yhat` is some collection of `DataFrame`s (e.g. Tuple, Vector).
"""
function idcol!(Yhat::Any, col::Symbol, val)
    for df in Yhat
        idcol!(df, col, val)
    end
end
