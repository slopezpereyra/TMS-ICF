using DataFrames
using CSV
using Statistics
using Plots
using StatsPlots
using LinearAlgebra
include("jl/math.jl")

df = CSV.read("df.csv", DataFrame)

function x_matrix(df::DataFrame, session_type::String)
    """Computes matrix k x n matrix X as defined in the research notes. X is 
    such that the row vectors r1, ..., rk are the paired responses of each of the 
    k experimental subjects."""

    filtered = filter([:ISI] => x -> x != -1, df)
    k = length(unique(df.Subject))
    X = Matrix{Float64}(undef, k, 192)

    for i in unique(df.Subject)
        subdf = filter([:Subject, :Session] => (x, y) -> x == i && y == session_type, filtered)
        if length(subdf.EMGPeakToPeak) != 192
            continue 
        end
        X = [X;transpose(subdf.EMGPeakToPeak)]
    end
    return X
end

function m_matrix(df::DataFrame, session_type::String, mprime::Bool=false)
    filtered = filter([:ISI, :Type] => (x, y) -> x == -1 && y == session_type, df)
    means::Vector{Float64} = []
    for i in unique(df.Subject)
        test_responses = filter([:Subject] => x -> x == i, filtered).EMGPeakToPeak
        if mprime 
            test_responses = test_responses.^(-1)
        end
        μ = mean(test_responses)
        push!(means, μ)
    end
    Diagonal(means)
end

X = x_matrix(df, "SWD")
M = m_matrix(df, "SWD")

M * X
