
function weighted_avg(a, penalties)
    sum(a .* penalties) / sum(penalties)
end

function p(d::Number, t::Vector, f=Φ)::Number
    """Computes weighted, pulse-specific relative
    amplitude given a paired-pulse d and a series
    of test pulses t. Weights are computed via some 
    penalty f."""
    w = broadcast(f, rzs(t))
    p = (d ./ t) .* w
    sum(p) / sum(w)
end

function p(d::Vector, t::Vector, f=Φ)::Vector
    h(x) = p(x, t, f)
    broadcast(h, d)
end

function p(d::SubArray, t::Vector, f=Φ)::Vector
    h(x) = p(x, t, f)
    broadcast(h, d)
end

function ℒ(θ::Vector, d::Vector, f)
    w = broadcast(f, d)
    sum(θ .* w) / sum(w)
end

# Penalty functions. See fourth chapter
# of paper in progress.

function Φ(x, α=3)
    α / (x^2 + α)
end

function Ψ(x, α=1)
    exp(-(x^2 / α))
end

function ∅(x)
    x
end

# Robust z-score

function rzs(x::Vector)::Vector
    """Given a vector x, returns a vector with the 
    robust z-scores of each element of x."""
    mad = median(abs.(x .- median(x)))
    0.675 * (x .- median(x)) / mad
end
