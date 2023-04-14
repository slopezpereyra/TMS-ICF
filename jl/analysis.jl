using DataFrames
using CSV
using Statistics
using Plots
using StatsPlots
include("jl/math.jl")

function sra(pulses, subject, session)
    """Computes standard relative amplitude given
    a subject's d vector on a given session."""
    # Using global variable
    t = get_test_pulses(subject, session)
    mean(pulses) / mean(t)
end

function get_test_pulses(subject, session)
    """Wrapper to avoid verbose expression"""
    filter([:ISI, :Subject, :Session] => (i, x, y) -> i == -1 && x == subject && y == session, df).EMGPeakToPeak
end

function subject_analysis(df)
    """Performs subject-level ICF/ICI analysis by computing the standard 
    relative amplitude on each session-subject subgroup per each ISI."""
    df = filter(:ISI => x -> x != 0 && x != -1, df)
    gdf = groupby(df, [:ISI, :Subject, :Session, :Label])
    an = combine(gdf, [:EMGPeakToPeak, :Subject, :Session] => (x, y, z) -> sra(x, y[1], z[1]))
    rename!(an, 5 => :SRA)
    sort(an, [:Subject, :Session, :ISI])
end

function group_analysis(an, col=:SRA)
    """Computes mean relative amplitude per ISI on each subject group."""
    an = groupby(an, [:Label, :ISI])
    an = combine(an, [col] => x -> mean(x))
    rename!(an, 3 => :SRA)
end

df = CSV.read("LnDF.csv", DataFrame)
display(df)    

s_df = subject_analysis(df)
display(s_df)
g_df = group_analysis(s_df)
display(df)


set_ρ(df)
set_ρ_ln(df)
set_ρ₂(df)
set_ρ₂_ln(df)
df = compute_inverse_variance_weights(df)
set_weighted_ρ()
set_weighted_ρ₂()

display(df)

CSV.write("full_df.csv", df)
