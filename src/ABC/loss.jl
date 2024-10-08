# loss.jl
# functions to compute the loss for two AbstractDatasets
# symmetric bounded loss is used by default
# observations should be configured using a yml config file


@enum SetPenalty proportional infinite

check_if_valid(x::AbstractVector) = @. isfinite(x) & !ismissing(x)

robust_log_transform(x::Real) = begin
    if x>= 0
        return (log(x + 1) + 1e-10)
    else
        return Inf
    end
end
robust_log_backtransform(x::Real) = exp(x - 1e10) - 1

# computes the kez
@inline length_based_weights(a::AbstractVector, b::AbstractVector) = Weights(
    vcat(
        fill(length(a), sum(@. !ismissing(a))), 
        fill(length(b), sum(@. !ismissing(b)))
        ))

# computes the weighted mean over two vectors, taking length differences into account
length_weighted_mean(a::AbstractVector, b::AbstractVector) = mean(vcat(a, b), length_based_weights(a, b))

"""
    symmetric_bounded_loss_with_penalty(
        predicted::AbstractVector,
        observed::AbstractVector;
        weights::Union{AbstractVector,Real} = 1,
        penalty_for_nans = proportional
        )

Symmetric bounded loss function  with penalty.  
This implementation of the Symmetric bounded loss function follows Marques et al. (2019), 
but applies a penalty to the loss if predictions contain non-finite values (including NaN and missing). 

The penalty is by default proportional to the fraction of predicted non-finite values. 
It can also be set to infinite, effectively rejecting the prediction altogether if there are any non-finite values in the prediction.

Predicted non-finite values can occur if we are exploring extreme parameter values and certain data points cannot be inferred from the simulation output. 
For example, if we fit a model to observed age at sexual maturity, it is possible that for some parameter samples maturity is never reached.

This function is used internally by the ABC fitting routine as inner loss function, and typically does not have to be called from the script or notebook. 
It can be helpful for diagnostic purposes, or as template for other loss functions.

args

- predicted: Vector of predicted values
- observed: Vector of observed values

kwargs

- weight: A single weight value
- penalty_for_nans: how to calculate the penalty for non-finited predicted values. `proportional` or `infinite`.
- transform: transformation applied to predicted and observed values. 
---

## Example

```Julia
using DEBBase.DEBODE, DEBBase.ABC, DataFrames

sim = DEBODE.simulator(Params())
pseudo_data = (7 .* (1 .- exp.(-0.2 .* sim.t))).^3

ABC.symmetric_bounded_loss_with_penalty(sim.S, pseudo_data)
```


"""
function symmetric_bounded_loss_with_penalty(
    predicted::AbstractVector,
    observed::AbstractVector;
    weight::Real = 1,
    penalty_for_nans = proportional,
    transform = x -> x + 1e-10
    )

    @assert length(predicted) == length(observed) "Got predicted and observed values of differing length: $(length(predicted)), $(length(observed))."

    valid_idcs_pred = check_if_valid(predicted)
    valid_idcs = valid_idcs_pred .& check_if_valid(observed)

    # this applies a penalty for missing values / NaNs in the predicted values
    penalty_factor = begin 
        if penalty_for_nans == proportional
            length(predicted)/sum(valid_idcs_pred)
        else
            Inf
        end
    end

    # filter vectors to only include non-missing and finite values
    obs_filt = @. transform(observed[valid_idcs])
    pred_filt = @. transform(predicted[valid_idcs])

    if isempty(pred_filt)
        return missing
    end

    let n = length(obs_filt), d = obs_filt, p = pred_filt
        pᵢ = (sum(p)/n)^2
        dᵢ = (sum(d)/n)^2
        loss = sum(@. (weight/n) * (( (d - p)^2 ) / ( pᵢ + dᵢ)))

        return loss * penalty_factor
    end
end


"""
    compute_time_resolved_loss(
        name::AbstractString,
        predicted::AbstractDataset, 
        observed::AbstractDataset; 
        inner_loss_function = symmetric_bounded_loss_with_penalty,
        return_all = false
        )::Union{Real,Missing}

Compute loss for a single time-resolved dataset entry.

args: 

- name: Name of the data entry in observed.time_resolved and predicted.time_resolved
- predicted: Predicted dataset
- observed: Observed dataset

kwargs:

- return_all: returns the losses for all response variables if true. Otherwise, return averaged loss. This option only works if directly calling this function, not if it is nested into compute_loss() or similar.
"""
function compute_time_resolved_loss(
    name::AbstractString,
    predicted::AbstractDataset, 
    observed::AbstractDataset; 
    weight = 1,
    inner_loss_function = symmetric_bounded_loss_with_penalty,
    return_all = false
    )

    observed_df = observed.time_resolved[name]
    predicted_df = predicted.time_resolved[name]
    time_var = observed.time_vars[name]
    grouping_vars = observed.grouping_vars["time_resolved"][name]
    response_vars = observed.response_vars["time_resolved"][name]

    # NOTE: this part causes most of the memory allocations in compute_loss, due to typeinf_ext_toplevel
    loss_per_response_var = observed_df |> 
    x -> leftjoin( # join predicted with observed 
        x, predicted_df, 
        on = vcat(time_var, grouping_vars), 
        makeunique = true, 
        renamecols = "_observed" => "_predicted"
        ) |> 
        drop_missing |>
        x -> groupby(x, grouping_vars) |> # apply desired grouping
        x -> combine(x) do df # calculate loss for each response variable

            losses = DataFrame()

            for response in response_vars
                col_observed = String(response)*"_observed" |> Symbol
                col_predicted = String(response)*"_predicted" |> Symbol
                loss = inner_loss_function(
                    df[:,col_predicted],
                    df[:,col_observed],
                    weight = weight
                    )

                append!(losses, DataFrame(response = response, loss = loss))
            end # for response

            return losses
        end # do df

    if return_all
        return loss_per_response_var
    end

    if isempty(skipmissing(loss_per_response_var.loss))
        return Inf
    end
    
    # return average and valid length of the data vector
    return mean(skipmissing(loss_per_response_var.loss))
end

"""
    compute_time_resolved_losses(
        predicted::AbstractDataset,
        observed::AbstractDataset;
        inner_loss_function = symmetric_bounded_loss_with_penalty
        )::Tuple

Computes loss values for all time-resolved data entries.
"""
function compute_time_resolved_loss(
    predicted::AbstractDataset,
    observed::AbstractDataset;
    weight = 1,
    return_all = false,
    inner_loss_function = symmetric_bounded_loss_with_penalty
    )

    datanames = keys(observed.time_resolved)
    time_resolved_losses = Union{Missing,Float64}[]

    for (i,name) in enumerate(datanames)
        loss = compute_time_resolved_loss(
            name, predicted, 
            observed; 
            inner_loss_function = inner_loss_function, 
            weight = observed.weights["time_resolved"][name]
            )
        push!(time_resolved_losses, loss)
    end

    if return_all
        return time_resolved_losses
    end

    return mean(skipmissing(time_resolved_losses))
end

# for scalar loss, we have three methods:
#   - one including the name of the dataset entry, dispatching to one of the other two
#   - one for data stored in dataframes
#   - one for data stored in dicts

# compute scalar loss for tabular data

"""
    compute_scalar_loss(
        predicted::AbstractDataFrame,
        observed::AbstractDataFrame,
        response_vars, 
        grouping_vars;
        inner_loss_function = symmetric_bounded_loss_with_penalty,
        return_all = false
        )

Computes loss values for scalar data in tabular format.
"""
function compute_scalar_loss(
    predicted::AbstractDataFrame,
    observed::AbstractDataFrame,
    response_vars, 
    grouping_vars;
    weight = 1,
    inner_loss_function = symmetric_bounded_loss_with_penalty,
    return_all = false
    )
    
    loss_per_response_var = observed |> 
    x -> begin # join predicted with observed
            if isempty(grouping_vars) 
                hcat(
                    rename(x, names(x) .* "_observed"),
                    rename(predicted, names(predicted) .* "_predicted")
                )
            else
                leftjoin(
                    x, predicted, 
                    on = grouping_vars,
                    makeunique = true, renamecols = "_observed" => "_predicted"
                    )
            end
        end |> 
        drop_missing |>
        x -> groupby(x, grouping_vars) |> # apply grouping
        x -> combine(x) do df # compute loss for each response var
            losses = DataFrame()

            for response in response_vars

                col_observed = String(response)*"_observed" |> Symbol
                col_predicted = String(response)*"_predicted" |> Symbol
                
                loss = inner_loss_function(
                    df[:,col_observed], 
                    df[:,col_predicted],
                    weight = weight
                    )
                append!(losses, DataFrame(
                    response = response, 
                    loss = loss
                    ), promote = true)
            end # for response

            return losses
        end

    if return_all
        return loss_per_response_var, nothing
    end

    if isempty(skipmissing(loss_per_response_var.loss))
        return Inf
    end

    return mean(skipmissing(loss_per_response_var.loss))
end

"""
    compute_scalar_loss(
        predicted::OrderedDict,
        observed::OrderedDict;
        inner_loss_function = symmetric_bounded_loss_with_penalty,
        return_all = false
        )

Computes loss values for scalar data in dicts (originates from yml files).
"""
function compute_scalar_loss(
    predicted::OrderedDict,
    observed::OrderedDict,
    response_vars,
    grouping_vars;
    weight = 1,
    inner_loss_function = symmetric_bounded_loss_with_penalty,
    return_all = false
    )

    scalar_losses = []

    for key in keys(observed)
        loss = inner_loss_function(
            [predicted[key]["value"]],
            [observed[key]["value"]],
            weight = weight
            ) 
        push!(scalar_losses, loss)
    end

    if return_all
        return Dict(zip(keys(observed), scalar_losses))
    end

    return mean(skipmissing(scalar_losses))
end


function compute_scalar_loss(
    predicted::AbstractDataset,
    observed::AbstractDataset;
    inner_loss_function = symmetric_bounded_loss_with_penalty
    )

    scalar_losses = []

    for name in keys(observed.scalar)
        push!(
            scalar_losses, 
            compute_scalar_loss(
                predicted.scalar[name],
                observed.scalar[name],
                observed.response_vars["scalar"][name],
                observed.grouping_vars["scalar"][name];
                weight = observed.weights["scalar"][name],
                inner_loss_function = inner_loss_function
            )
        )
    end

    return mean(skipmissing(scalar_losses))
end

"""
    compute_loss(
        predicted::AbstractDataset,
        observed::AbstractDataset;
        inner_loss_function = symmetric_bounded_loss_with_penalty
        )

Computes the loss from two datasets. 
Requires configuration of the calibration data, e.g. through a config file 
(cf test/config/data_config_example.yml).

Used internally by `DEBBase.ABC`, but can be used for customized calibration routines.
"""
function compute_loss(
    predicted::AbstractDataset,
    observed::AbstractDataset;
    return_all = false,
    inner_loss_function = symmetric_bounded_loss_with_penalty
    )

    time_resolved_loss = compute_time_resolved_loss(
        predicted, 
        observed; 
        inner_loss_function = inner_loss_function
        )
    
    scalar_loss = compute_scalar_loss(
        predicted,
        observed;
        inner_loss_function = inner_loss_function
    )
    
    if return_all
        return time_resolved_loss, scalar_loss
    end

    total_loss = mean([time_resolved_loss, scalar_loss])

    # if there are still missing values, this leads to infite distance, 
    # which in turn leads to rejection of the sample in the SMC algorithm
    if ismissing(total_loss)
        return Inf
    else
        return total_loss
    end 
end