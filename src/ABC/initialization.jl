"""
Extract fields from a parameter structure.
"""
function extractfields(p::Any)
    return [getproperty(p, f) for f in fieldnames(typeof(p))]
end


"""
    deftruncnorm(mu, CV; l = 0, u = Inf)

Define a truncated normal distribution based on the mean and CV.
"""
function deftruncnorm(mu, CV; l = 0, u = Inf) 
    return Truncated(Normal(mu, abs(CV * mu)), l, u)
end

"""
    deflognorm(modus, sigma)

Define a truncated normal distribution based on the mode and σ. 
The parameter μ is then calculated as ``\\mu = ln(mode) + \\sigma^2``.
"""
function deflognorm(modus, sigma)
    mu = log(modus) + sigma^2
    return LogNormal(mu, sigma)
end