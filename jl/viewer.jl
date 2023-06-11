
using DataFrames
using CSV
using Statistics
using Plots
using StatsPlots

df2 = CSV.read("/data/SWIP_009_S1.txt", delim='\')

test = filter(row -> row.Subject == 9 && row.Session == 1, df)
display(test)

plot(test.EMGPeakToPeak)

