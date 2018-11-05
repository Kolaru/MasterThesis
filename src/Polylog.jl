module Polylog

export polylog, polylog_over_z, dpolylog_over_z

using PyCall

const mpmath = PyNULL()
const polylog_obj = PyNULL()

"""
The way to import `polylog` from Python `mpmath` library is convoluted, because we
want to use the fast floating point implementation with the `fp` (floating point)
context.
"""
function __init__()
    copy!(mpmath, pyimport("mpmath"))
    copy!(polylog_obj, mpmath[:functions][:functions][:SpecialFunctions][:defined_functions]["polylog"][1])
end

"""
    polylog(s, z)

Polylogarithm function.
"""
polylog(s::Real, z::Real)::Float64 = abs(polylog_obj(mpmath[:fp], s, z))

function polylog_over_z(s::Real, z::Real)
    z < 1e-8 && return 1. + z/2^s + z^2/3^s
    return polylog(s, z)/z
end

function dpolylog_over_z(s::Real, z::Real)
    if z == 0
        return 1/(2^s)
    else
        return (polylog_over_z(s-1, z) - polylog_over_z(s, z))/z
    end
end

end
