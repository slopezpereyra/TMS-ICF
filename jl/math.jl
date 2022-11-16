function 𝒫(𝓭::Number, 𝓽::Vector, f)::Number
    """Computes weighted, pulse-specific relative
    amplitude given a paired-pulse d and a series
    of test pulses t. Weights are computed via some 
    penalty f."""
    𝔀 = broadcast(f, rzs(𝓽))
    p = (𝓭 ./ 𝓽) .* 𝔀
    sum(p) / sum(𝔀)
end

function ℒ(θ::Vector, d::Vector, f)
    w = broadcast(f, d)
    sum(θ .* w) / sum(w)
end

function Φ(x, α=1)
    α / (x^2 + α)
end

function rzs(x::Vector)::Vector
    mad = median(abs.(x .- median(x)))
    0.675 * (x .- median(x)) / mad

end
