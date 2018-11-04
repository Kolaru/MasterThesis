__precompile__()

module GeneratingFunctions

using Graphs
using MonotonicExtension

import LerchPhi: lerchphi
import Polylog: polylog, polylog_over_z, dpolylog_over_z

using IntervalArithmetic
using PyCall

import SpecialFunctions: zeta

export g0, dg0, g1, dg1
export polylog, polylog_over_z

const ZETAS = Dict{Any, Float64}()
const MAX_EXP = 100

"""
    zeta_storing(α)

Return the same result as the `zeta(α)` function but store every value computed.
More efficient for the cases where the value of `zeta(α)` is requested multiple
times for the same `\alpha`.
"""
function zeta_storing(args...)
    if !haskey(ZETAS, args)
        ZETAS[args] = zeta(args...)
    end
    return ZETAS[args]
end

# Extend polylog family for interval arithmetic
@monotone polylog Domain(2..MAX_EXP, 0..1)
@monotone polylog_over_z Domain(2..MAX_EXP, 0..1)
@monotone dpolylog_over_z Domain(2..MAX_EXP, 0..1)

@monotone zeta Domain(2..MAX_EXP)
@monotone zeta Domain(2..MAX_EXP, 0..1000)

# Common interface for all network types
g0(::Type{ErdosRenyiGraph}, z, c) = exp(c*(z-1))
dg0(::Type{ErdosRenyiGraph}, z, c) = c * g0(ErdosRenyiGraph, z, c)
g1(::Type{ErdosRenyiGraph}, z, c) = g0(ErdosRenyiGraph, z, c)
dg1(::Type{ErdosRenyiGraph}, z, c) = dg0(ErdosRenyiGraph, z, c)

# abs(.) is to avoid imaginary number when the functions are called outside of
# their allowed range (typically done by nlsolve)
g0(::Type{ScaleFreeGraph}, z, α) = abs(polylog(α, z)/zeta(α))
dg0(::Type{ScaleFreeGraph}, z, α) = abs(polylog_over_z(α-1, z)/zeta(α))
g1(::Type{ScaleFreeGraph}, z, α) = abs(polylog_over_z(α-1, z)/zeta(α-1))
dg1(::Type{ScaleFreeGraph}, z, α) = abs(dpolylog_over_z(α-1, z)/zeta(α-1))

# Using @monotone on g0 and g1 directly would have been better than on
# the zeta and polylog family. This works, so let it be for now.
function g0(::Type{ScaleFreeGraph}, z::Interval, α)
    res = polylog(α, z)/zeta(α)
    return Interval(max(res.lo, 0), min(res.hi, 1))
end

function g1(::Type{ScaleFreeGraph}, z::Interval, α)
    res = polylog_over_z(α-1, z)/zeta(α-1)
    return Interval(max(res.lo, 0), min(res.hi, 1))
end

ssf_g1(z, s, a) = (lerchphi(z, s-1, a+1) - a*lerchphi(z, s, a+1))/(zeta_storing(s-1, a+1) - a*zeta_storing(s, a+1))

# @monotone ssf_g1 IntervalBox(0..1, 2..MAX_EXP, 0..1000) 1e-11 0..1

g0(::Type{SaturatedScaleFreeGraph}, z, s, a) = z*lerchphi(z, s, a+1)/zeta(s, a+1)
dg0(::Type{SaturatedScaleFreeGraph}, z, s, a) = (lerchphi(z, s-1, a+1) - a*lerchphi(z, s, a+1))/zeta(s, a+1)
g1(::Type{SaturatedScaleFreeGraph}, z, s, a) = ssf_g1(z, s, a)
dg1(::Type{SaturatedScaleFreeGraph}, z, s, a) = error("not implemented")

g0(::Type{GeometricGraph}, z, c) = 1/(c*(1/z - 1) + 1)
dg0(::Type{GeometricGraph}, z, c) = c/(c*(1 - z) + z)^2
g1(::Type{GeometricGraph}, z, c) = 1/(c*(1 - z) + z)^2
dg1(::Type{GeometricGraph}, z, c) = 2 * (c - 1) * 1/(c - (c - 1)*z)^3

# Special definition to avoid the problem X/X != 1
function g1(::Type{GeometricGraph}, z::Interval, c)
    A = c*(1 - z)
    lim1 = A.lo + z.hi
    lim2 = A.hi + z.lo
    a = min(lim1, lim2)
    b = max(lim1, lim2)
    return 1/(a..b)^2
end

for func in (:g0, :dg0, :g1, :dg1)
    # Definitions of the closures
    @eval ($func)(::Type{T}) where {T <: GraphType} = (args...) -> ($func)(T, args...)
end

include("GeneratingFunctions/multilayer.jl")

end
