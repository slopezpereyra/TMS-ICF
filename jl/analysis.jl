using DataFrames
using CSV
using Statistics
using Plots
include("math.jl")

df = DataFrame(CSV.File("TMS-ICF/df.csv"))
df = select!(df, Not(:Column1))
df = sort(df, [:Subject, :Session, :ISI])

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

function set_pulse_ara(f, name::Symbol)
    # Not working
    gdf = groupby(df, [:ISI, :Subject, :Session])
    res = combine(gdf, [:EMGPeakToPeak, :Subject, :Session] => (x, y, z) -> p(x, get_test_pulses(y[1], z[1]), f))
    res = sort(res, [:Subject, :Session, :ISI])
    df[!, name] = res[:, end]
end

function set_pulse_rra()
    gdf = groupby(df, [:ISI, :Subject, :Session])
    res = combine(gdf, [:EMGPeakToPeak, :Subject, :Session] => (x, y, z) -> rawp(x, get_test_pulses(y[1], z[1])))
    res = sort(res, [:Subject, :Session, :ISI])
    df[!, :_RRA] = res[:, end]
end

function subject_analysis(df)
    df = filter(:ISI => x -> x != 0 && x != -1, df)
    df = groupby(df, [:ISI, :Subject, :Session, :Label])
    an = combine(df, [:EMGPeakToPeak, :Subject, :Session] => (x, y, z) -> sra(x, y[1], z[1]))
    rename!(an, 5 => :SRA)
    sort(an, [:Subject, :Session, :ISI])
end

function group_analysis(an)
    an = groupby(an, [:Label, :ISI])
    an = combine(an, [:SRA] => x -> mean(x))
    rename!(an, 3 => :SRA)
end

df

function set_pulse_specific_ras()
    set_pulse_rra()
    set_pulse_ara(Φ, :QRA)
    set_pulse_ara(Ψ, :GaussRA)
end

set_pulse_specific_ras()
an = subject_analysis(df)
gan = group_analysis(an)