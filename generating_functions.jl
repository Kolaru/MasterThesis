include("polylog.jl")

const ZETAS = Dict{Float64, Float64}()

"""
    zeta_storing(α)

Return the same result as the `zeta(α)` function but store every value computed.
More efficient for the cases where the value of `zeta(α)` is requested multiple
times for the same `\alpha`.
"""
function zeta_storing(α::Number)
    if !haskey(ZETAS, α)
        ZETAS[α] = zeta(α)
    end
    return ZETAS[α]
end

function polylog_over_z(s::Real, z::Real)
    z == 0 && return zero(z)
    return polylog(s, z)/z
end

function dpolylog_over_z(s::Real, z::Real)
    if z == 0
        return 1/(2^s)
    else
        return (polylog_over_z(s-1, z) - polylog_over_z(s, z))/z
    end
end

# Extend polylog family for interval arithmetic
@extend_monotonic polylog(:down, :up)
@extend_monotonic polylog_over_z(:down, :up)
@extend_monotonic dpolylog_over_z(:down, :up)

# Common interface for all network types
g0(::Type{ErdosRenyiGraph}, z, c) = exp(c*(z-1))
dg0(::Type{ErdosRenyiGraph}, z, c) = c * g0(ErdosRenyiGraph, z, c)
g1(::Type{ErdosRenyiGraph}, z, c) = g0(ErdosRenyiGraph, z, c)
dg1(::Type{ErdosRenyiGraph}, z, c) = dg0(ErdosRenyiGraph, z, c)

g0(::Type{ScaleFreeGraph}, z, α) = polylog(α, z)/zeta_storing(α)
dg0(::Type{ScaleFreeGraph}, z, α) = polylog_over_z(α-1, z)/zeta_storing(α)
g1(::Type{ScaleFreeGraph}, z, α) = polylog_over_z(α-1, z)/zeta_storing(α-1)
dg1(::Type{ScaleFreeGraph}, z, α) = dpolylog_over_z(α-1, z)/zeta_storing(α-1)

g0(::Type{GeometricGraph}, z, c) = 1/c * z/(1 - (1 - 1/c)*z)
dg0(::Type{GeometricGraph}, z, c) = 1/c * 1/(1 - (1 - 1/c)*z)^2
g1(::Type{GeometricGraph}, z, c) = dg0(GeometricGraph, z, c)/c
dg1(::Type{GeometricGraph}, z, c) = error("Not implemented yet")

# Definitions of the closures
for func in (:g0, :dg0, :g1, :dg1)
    @eval ($func)(::Type{T}) where {T <: GraphGenerator} = (z, c) -> ($func)(T, z, c)
end
