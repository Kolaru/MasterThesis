module Specials

using IntervalArithmetic
using SpecialFunctions

import SpecialFunctions.zeta

export polylog, polylog_over_z, zeta

zeta(s::Real, z::Interval{T}) where {T <: Real} = widen(zeta(s, z.hi)..zeta(s, z.lo))

function zeta(s::Real, z::Complex{Interval{T}}) where {T <: Real}
    println(z)
    error()
end


function polylog_old(s::Real, z::Real)
    ints = round(Int, s)
    if isapprox(s, ints, atol=1e-6)
        return polylog(ints, z)
    end

    twopi = 2π
    x = log(z)/(im*twopi)
    ss = 1 - s
    ip = im^(1-s)
    return real(gamma(ss)/twopi^ss * (ip*zeta(ss, 1+x) + conj(ip)*zeta(ss, -x)))
end


function polylog_over_z_naive(s::Real, z::Real)
    if z > 1
        error("Real polylogarithm does not converge for `z > 1`.")
    end
    tol = 1e-10
    sum = 0
    zk = 1
    k = 1
    term = 1

    while term > tol
        term = zk / k^s
        sum += term
        zk *= z
        k += 1
    end

    return sum
end

const LOGTOL = log(eps(1.0))

function polylog_over_z(s::Real, z::Real)
    # println("Polylog real")
    # @show s
    # @show z

    if z == 1
        return zeta(s)
    end

    if z > 1
        error("Real polylogarithm does not converge for `z > 1`.")
    end

    K = ceil(LOGTOL/log(z))
    sum = 0.0

    # for k in K:-1:1
    #     sum += z^(k-1)/k^s
    # end

    while K >= 1
        sum += z^(K-1)/K^s
        K -= 1
    end
    return sum
end


# function polylog_over_z(s::Real, z::Interval{T}) where {T <: Real}
#     println("Polylog interval")
#     @show s
#     @show z
#     if z.hi > 1
#         error("Real polylogarithm does not converge for `z > 1`.")
#     end
#     K = ceil(LOGTOL/log(z.lo))
#     sum = 0.0
#
#     for k in K:-1:1
#         sum += z^(k-1)/k^s
#     end
#     return sum
# end



polylog(s, z) = z*polylog_over_z(s, z)


end
