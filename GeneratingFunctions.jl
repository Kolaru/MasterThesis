module GeneratingFunctions

using SpecialFunctions
using Specials

export Geometric, Poisson, Powerlaw


function poisson_g0(z, c)
    return exp(c*(z - 1))
end

function poisson_dg0(z, c)
    return c * exp(c*(z - 1))
end

function powerlaw_g0(z, a)
    return polylog(a, z)/zeta(a)
end

function powerlaw_g1(z, a)
    return polylog_over_z(a-1, z)/zeta(a-1)
end

function powerlaw_dg1(z, a)
    return (polylog_over_z(a-2, z)/zeta(a-1) - powerlaw_g1(z, a))/z
end

function geometric_g0(z, c)
    return 1/c * z/(1 - (1 - 1/c)*z)
end

function geometric_dg0(z, c)
    return 1/c * 1/(1 - (1 - 1/c)*z)^2
end

function geometric_g1(z, c)
    return geometric_dg0(z, c)/c
end


abstract type DegreeDistribution{G0, G1, DG1} end


struct Poisson{G0, G1, DG1} <: DegreeDistribution{G0, G1, DG1}
    g0::G0
    g1::G1
    dg1::DG1
end

function Poisson(c::Real)
    g0 = z -> poisson_g0(z, c)
    g1 = z -> poisson_g0(z, c)
    dg1 = z -> poisson_dg0(z, c)
    Poisson(g0, g1, dg1)
end


struct Powerlaw{G0, G1, DG1} <: DegreeDistribution{G0, G1, DG1}
    g0::G0
    g1::G1
    dg1::DG1
end

function Powerlaw(c::Real)
    g0 = z -> powerlaw_g0(z, c)
    g1 = z -> powerlaw_g1(z, c)
    dg1 = z -> powerlaw_dg1(z, c)
    Powerlaw(g0, g1, dg1)
end


struct Geometric{G0, G1, DG1} <: DegreeDistribution{G0, G1, DG1}
    g0::G0
    g1::G1
    dg1::DG1
end

function Geometric(c::Real)
    g0 = z -> geometric_g0(z, c)
    g1 = z -> geometric_g1(z, c)
    Geometric(g0, g1, dg1)
end

end
