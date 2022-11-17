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


function rawp(d::Vector, t::Vector)::Vector
    d ./ mean(t)
end

function rawp(d::SubArray, t::Vector)::Vector
    d ./ mean(t)
end

function ℒ(θ::Vector, d::Vector, f)
    w = broadcast(f, d)
    sum(θ .* w) / sum(w)
end

function Φ(x, α=1)
    α / (x^2 + α)
end

function Ψ(x, α=1)
    exp(-(x^2 / α))
end

function ∅(x)
    x
end
function rzs(x::Vector)::Vector
    mad = median(abs.(x .- median(x)))
    0.675 * (x .- median(x)) / mad

end
