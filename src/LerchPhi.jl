module LerchPhi

using Cubature

export lerchphi

const twopi = 2*pi
const QUAD_RELTOL = 1e-12

second_term(z, s, a, t) = z^t/(a + t)^s
third_term(alogz, s, a, t) = (sin(s*atan(t) - t*alogz))/((1 + t^2)^(s/2)*(exp(twopi*t*a) - 1))

zero_to_infinity(u) = u/(1-u)

function lerchphi_raw(z, s, a)
    z == 0 && return a^(-s)

    res = 1/(2*a^s)
    sterm(t) = second_term(z, s, a, zero_to_infinity(t)) / (1-t)^2
    (sval, serr) = hquadrature(sterm, 0, 1 ; reltol=QUAD_RELTOL)
    res += sval
    alogz = a*log(z)
    tterm(t) = third_term(alogz, s, a, zero_to_infinity(t)) / (1-t)^2
    (tval, terr) = hquadrature(tterm, 0, 1 ; reltol=QUAD_RELTOL)
    return res + 2/(a^(s-1))*tval
end

const LERCHPHIS = Dict()
function lerchphi_storing(args...)
    if !haskey(LERCHPHIS, args)
        LERCHPHIS[args] = lerchphi_raw(args...)
    end
    return LERCHPHIS[args]
end


lerchphi = lerchphi_storing

end
