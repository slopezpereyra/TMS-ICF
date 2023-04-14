using DataFrames
using CSV
using Statistics
using Plots
using StatsPlots
using Distributions

include("jl/math.jl")

# TEST AND NO-TEST DISTRIBUTIONS 

qqplot(Exponential(248), df.EMGPeakToPeak)
mean(df.EMGPeakToPeak)

no_tests = filter(:ISI => x -> x != -1, df).EMGPeakToPeak
tests = filter(:ISI => x -> x == -1, df).EMGPeakToPeak
mean(no_tests)

qqplot(Exponential(260), no_tests)

density(no_tests, linewidth=2, color=:blue, label="Paired pulse distribution")
density!(tests, linewidth=2, color=:red, label = "Test pulse distribution")
xlabel!("EMG Peak to peak")
ylabel!("Probability density")

# GROUP - TYPE DISTRIBUTIONS 

function get_group(group, type, tests)
    data = filter([:Group, :Type] => (x, y) -> x == group && y == type, df)
    if tests 
        data = filter([:ISI] => x -> x == -1, data).EMGPeakToPeak
    else
        data = filter([:ISI] => x -> x != -1, data).EMGPeakToPeak
    end
    return data
end 

HCBL = get_group(1, "BL", false)
HCSWD = get_group(1, "SWD", false)
MDDBL = get_group(2, "BL", false)
MDDSWD = get_group(2, "SWD", false)


density(HCBL, linewidth=2, color=:blue, label="HC BL")
density!(HCSWD, linewidth=2, color=:green, label = "HC SWD")
density!(MDDBL, linewidth=2, color=:red, label="MDD BL")
density!(MDDSWD, linewidth=2, color=:purple, label = "MDD SWD")
title!("Paired pulse distribution across subject groups")
xlabel!("EMG Peak to peak")
ylabel!("Probability density")

median(HCBL)
median(HCSWD)
median(MDDBL)
median(MDDSWD)

