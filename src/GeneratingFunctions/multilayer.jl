export SingleParameterGraph
export ψ, ψ!, dψ

struct SingleParameterGraph{T} <: GraphType
    L::Int  # Number of layers
end

"""
    Fixpoint function `ψ`.

Order of arguments `λ` and `z` are inverted compare to what is done in the
report, to match the definition of `g0` and `g1` in `GeneratingFunctions.jl`.
"""
function ψ(sp::SingleParameterGraph{GT}, z, λ) where {GT <: GraphType}
    return 1. - (1. - g1(GT, z, λ))*(1 - g0(GT, z, λ))^(sp.L - 1)
end

function ψ!(layers::Vector{DataType}, z, λ, res)
    L = length(layers)
    p0 = [1 - g0(layers[k], z[k], λ[k]) for k in 1:L]
    p1 = [1 - g1(layers[k], z[k], λ[k]) for k in 1:L]

    for k in 1:L
        prod = 1.
        for j in 1:L
            if k == j
                prod *= p1[j]
            else
                prod *= p0[j]
            end
        end
        res[k] = 1 - prod
    end
    return res
end

function ψ(layers::Vector{DataType}, z::IntervalBox{N, T}, λ) where {N, T}
    return IntervalBox(ψ!(layers, z, λ, Vector{Interval{T}}(undef, N)))
end

function ψ(layers::Vector{DataType}, z::T, λ) where {T}
    return ψ!(layers, z, λ, Vector{T}(undef, length(layers)))
end

ψ(layers::Vector{DataType}) = (z, λ) -> ψ(layers, z, λ)

function dψ(sp::SingleParameterGraph{GT}, z, λ) where {GT <: GraphType}
    p0 = 1 - g0(GT, z, λ)
    return ( dg1(GT, z, λ)*p0 + (sp.L - 1)*(1 - g1(GT, z, λ))*dg0(GT, z, λ) ) * p0^(sp.L-2)
end
