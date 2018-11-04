export SingleParameterGraph
export ψ, ψ!, dψ, dψ!

struct SingleParameterGraph{T} <: GraphType
    L::Int  # Number of layers
end

"""
    Fixpoint function `ψ`.

Order of arguments `λ` and `z` are inverted compare to what is done in the
report, to match the definition of `g0` and `g1` in `GeneratingFunctions.jl`.

It is assumed that each layer is dependent on exactly one parameter.
"""
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
        res[k] = 1. - prod
    end
    return res
end

function ψ(layers::Vector{DataType}, z::IntervalBox{N, T}, λ) where {N, T}
    res = Vector{Interval{T}}(undef, N)
    ψ!(layers, z, λ, res)
    return IntervalBox(res)
end

function ψ(layers::Vector{DataType}, z::Vector{T}, λ) where {T}
    res = Vector{T}(undef, length(layers))
    ψ!(layers, z, λ, res)
    return res
end

function ψ(sp::SingleParameterGraph{GT}, z, λ) where {GT <: GraphType}
    return 1. - (1. - g1(GT, z, λ))*(1 - g0(GT, z, λ))^(sp.L - 1)
end

function dψ!(layers::Vector{DataType}, z, λ, J)
    L = length(layers)
    p0 = [1. - g0(layers[k], z[k], λ[k]) for k in 1:L]
    p1 = [1. - g1(layers[k], z[k], λ[k]) for k in 1:L]
    dp0 = [-dg0(layers[k], z[k], λ[k]) for k in 1:L]
    dp1 = [-dg1(layers[k], z[k], λ[k]) for k in 1:L]

    for i in 1:L, j in 1:L
        prod = 1.
        for k in 1:L
            if k == i
                if i == j
                    prod *= dp1[k]
                else
                    prod *= p1[k]
                end
            elseif k == j
                prod *= dp0[k]
            else
                prod *= p0[k]
            end
        end

        J[i, j] = -prod
    end
end

function dψ(layers::Vector{DataType}, z::IntervalBox{N, T}, λ) where {N, T}
    J = Array{Interval{T}, 2}(undef, N, N)
    dψ!(layers, z, λ, J)
    return J
end

function dψ(layers::Vector{DataType}, z::Vector{T}, λ) where {T}
    L = length(layers)
    J = Array{T, 2}(undef, L, L)
    dψ!(layers, z, λ, J)
    return J
end

function dψ(sp::SingleParameterGraph{GT}, z, λ) where {GT <: GraphType}
    p0 = 1. - g0(GT, z, λ)
    return ( dg1(GT, z, λ)*p0 + (sp.L - 1)*(1. - g1(GT, z, λ))*dg0(GT, z, λ) ) * p0^(sp.L-2)
end

for func in (:ψ, :dψ)
    @eval ($func)(layers::Vector{DataType}) = (z, λ) -> ψ(layers, z, λ)
end
