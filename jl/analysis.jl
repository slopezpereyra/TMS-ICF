using DataFrames
using CSV
using Statistics

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

function pulse_ra()
    # Not working
    gdf = groupby(df, [:ISI, :Subject, :Session])
    combine(gdf, [:EMGPeakToPeak, :Subject, :Session] => (x, sub, ses) -> 𝒫(x, get_test_pulses(sub, ses)))
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
