"""
Extract fields from a parameter structure.
"""
function extractfields(p::Any)
    return [getproperty(p, f) for f in fieldnames(typeof(p))]
end


"""
Define a truncated normal distribution based on the mean and CV.
"""
function deftruncnorm(mu, CV; l = 0, u = Inf)  # TODO: move deftruncnorm() to DEBABC
    return Truncated(Normal(mu, CV * mu), l, u)
end

function deflognorm(modus, sigma)
    mu = log(modus) + sigma^2
    return LogNormal(mu, sigma)
end