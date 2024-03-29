using CSV
using Statistics
using DataFrames

# Inverse variance functions

function single_point_variance(x, vec)
    """Compute single-point variance of x with respect to vector
    of origin. It must be the case that x ∈ vec"""
    filtered = filter(e -> e != x, vec)
    sum(abs.(x .- filtered)) / length(filtered) # Length filtered = n - 1 since x was removed.
end

function inv_var_weights(vec)
    """Compute inverse-variance weights of all elements
    in a given vector"""
    f(x) = single_point_variance(x, vec) .^ (-1)
    f.(vec)
end

function compute_inverse_variance_weights(df)
    """Wrapper function that adds the corresponding inverse-variance weight 
    to each pulse in the data."""
    t = groupby(df, [:ISI, :Label])
    df = transform(t, :EMGPeakToPeak .=> inv_var_weights => :InvVarWeights)
end # Helper functions

function get_pulses(subject, session, isi, colname=:EMGPeakToPeak)
    """Returns some column of the data filtered by subject, session, and isi. 
    Defaults to EMG Peak to peak."""
    filter([:ISI, :Subject, :Session] => (i, x, y) -> i == isi && x == subject && y == session, df)[!, colname]
end

function get_test_pulses(df, subject, session, colname=:EMGPeakToPeak)
    """Wrapper to avoid verbose expression"""
    filter([:ISI, :Subject, :Session] => (i, x, y) -> i == -1 && x == subject[1] && y == session[1], df)[!, colname]
end
# Pulse-specific relative amplitudes

function ρ(d::Vector, t::Vector)::Vector
    d ./ mean(t)
end

function ρ(d::SubArray, t::Vector)::Vector
    d ./ mean(t)
end

function ρ₂(d, t)::Vector
    f(x) = mean(x .* (t .^ -1))
    broadcast(f, d)
end

function weighted_ρ(d, t, weights)
    """Compute the relative amplitude of a given paired-pulse d
    with respect to the weighted mean of t"""
    wt = weighted_avg(t, weights)
    d / wt
end

function weighted_ρ₂(d, t, weights)
    """Computes second form of the weighted relative amplitude. 
    See notes."""
    f(x) = weighted_avg(x ./ t, weights)
    broadcast(f, d)
end

### Setters

function set_weighted_ρ()
    """Sets weghted relative amplitude feature in the database."""
    gdf = groupby(df, [:ISI, :Subject, :Session])
    res = transform(gdf, [:EMGPeakToPeak, :Subject, :Session] => (x, y, z) -> weighted_ρ(x, get_test_pulses(df, y, z), get_test_pulses(df, y, z, :InvVarWeights)))
    df[!, :WRA] = res[:, end]
end

function set_weighted_ρ₂()
    """Sets weighted relative amplitude feature in the database using the
    second form of the weighted relative amplitude."""
    gdf = groupby(df, [:ISI, :Subject, :Session])
    res = transform(gdf, [:EMGPeakToPeak, :Subject, :Session] => (x, y, z) -> weighted_ρ₂(x, get_test_pulses(df, y, z), get_test_pulses(df, y, z, :InvVarWeights)))
    df[!, :WRA2] = res[:, end]
end


function set_ρ(df)
    """Sets raw, pulse-specific relative amplitude feature in the database."""
    gdf = groupby(df, [:ISI, :Subject, :Session])
    f(x, y, z) = ρ(x, get_test_pulses(df, y, z))
    res = transform(gdf, [:EMGPeakToPeak, :Subject, :Session] => (x, y, z) -> ρ(x, get_test_pulses(df, y, z)))
    df[!, :RRA] = res[:, end]
end

function set_ρ_ln(df)
    """Sets raw, pulse-specific relative amplitude feature in the database."""
    gdf = groupby(df, [:ISI, :Subject, :Session])
    f(x, y, z) = ρ(log.(x), log.(get_test_pulses(df, y, z)))
    res = transform(gdf, [:EMGPeakToPeak, :Subject, :Session] => (x, y, z) -> f(x, y, z))
    df[!, :LogRRA] = res[:, end]
end


function set_ρ₂(df)
    """Sets raw, pulse-specific relative amplitude feature in the database."""
    gdf = groupby(df, [:ISI, :Subject, :Session])
    f(x, y, z) = ρ₂(x, get_test_pulses(df, y, z))
    res = transform(gdf, [:EMGPeakToPeak, :Subject, :Session] => (x, y, z) -> f(x, y, z))
    df[!, :RRA2] = res[:, end]
end

function set_ρ₂_ln(df)
    """Sets raw, pulse-specific relative amplitude feature in the database."""
    gdf = groupby(df, [:ISI, :Subject, :Session])
    f(x, y, z) = ρ₂(log.(x), log.(get_test_pulses(df, y, z)))
    res = transform(gdf, [:EMGPeakToPeak, :Subject, :Session] => (x, y, z) -> f(x, y, z) )
    df[!, :LogRRA2] = res[:, end]
end
