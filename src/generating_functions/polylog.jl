# Taken from https://github.com/JuliaMath/SpecialFunctions.jl/pull/30
# Author: Matthew Roughan
# GitHub account: https://github.com/mroughan
# TODO Give credit and send thanks message
# Replace all that by:
# mpmath.functions[:functions][:SpecialFunctions][:defined_functions]["polylog"][1]

include("harmonic.jl")

function polylog(s::Number, z::Number, accuracy::Real=1.0e-18)
    T = 0.5 # threshold at which we change algorithms
    if z ≈ 1.0
        if real(s) > 1
            return zeta(s)
        else
            return Inf
        end
    elseif z ≈ -1.0
        return -eta(s)
    elseif s ≈ 0.0
        return z ./ (1-z)
    elseif abs(z) <= T
        return polylog_direct(s, z, accuracy)
    elseif abs(z) >= 1/T && isinteger(s) && real(s) < 0
        # use inversion formula to calculate in terms of Li_n(1/z)
        # but note for negative integer s, it collapses to something small
        return -(-1.0)^s .*polylog_direct(s, 1/z, accuracy)
    elseif  abs(z) >= 1/T
        # use inversion formula to calculate in terms of Li_s(1/z)
        z = convert(Complex{Float64}, z)
        G = (2*pi*im)^s * zeta( 1-s, 0.5 + log(-z)/(2*pi*im) ) /  gamma(s)
        F = complex(-1.0)^s * polylog_direct(s, 1/z, accuracy)

        A = 2*pi*im*log(z)^(s-1)/(gamma(s))
        if ( isreal(z) && real(z)>=1 )
            Θ = 1
        else
            Θ = 0
        end
        # println("G = $G, F=$F, Θ=$Θ, A=$A")
        return ( G - F - Θ*A )
    else
        # power series around mu=0, for z = e^mu
        polylog_series_mu(s, z, accuracy)
    end
end


####################################################################
#### these are component functions and aren't exported at this point

# Dirichlet beta function, for testing results
#   https://en.wikipedia.org/wiki/Dirichlet_beta_function
#   but don't export as it would create confusion with Euler integral
function Dbeta(s::Number)
    β = 4.0^(-s) * ( zeta(s,0.25) - zeta(s,0.75) )
end

function polylog_zeta(s::Number, z::Number, accuracy=1.0e-18)
    # compute using the Hurwitz-zeta function identity
    twopi = 2π
    x = im * (log(complex(-z)) / twopi)
    ss = 1-s
    ip = im^ss
    return ( gamma(ss)/twopi^(ss) ) * (ip * zeta(ss, 0.5-x) + conj(ip) * zeta(ss, 0.5+x))
end


function polylog_direct(s::Number, z::Number, accuracy=1.0e-18)
    # calculate using direct definition
    if abs(z) > 1 || ( abs(z) ≈ 1  && real(s) <= 2)
        throw(DomainError())
    end
    if abs(z) > 1/2
        warn("Slow convergence for  |z| > 1/2")
    end
    total = 0.0
    L = ceil(-log10(accuracy)*log2(10)) # summation limit from Crandall, which is conservative, but based on real s
    a = z;
    for n=1:L
        total += a
        a *= z * ( n/(n+1.0) )^s
        # println("   total = $total")
        if abs(a)/abs(total) < 1.0e-30
            break
        end
    end
    return total
end

function polylog_series_mu(s::Number, z::Number, accuracy=1.0e-18)
    # calculate using power series around μ = log(z) = 0
    μ = log(convert(Complex{Float64}, z))
    # println("μ = $μ")
    if abs(μ) > 2*pi
        throw(DomainError())
    end
    L = Int(ceil(-log10(accuracy)*log2(10))) # revisit this limit
    if isinteger(s)
        n = Int(round(s))
        if n>1
            # Crandall's 1.4 for s integer
            total = μ^(n-1)*(harmonic(n-1) - log(-μ))/gamma(n)
            # println("   μ=$μ, total = $total")
            tmp = 1
            for m=0:L
                if n - m != 1
                    # total += μ^m * zeta(n - m) / gamma(m+1)
                    total += tmp * zeta(n - m)
                end
                # println("   m=$m, total = $total, tmp=$tmp, ctmp=$(μ^m /gamma(m+1))")
                tmp *= μ/(m+1)
                if abs(tmp)/abs(total) < 1.0e-30
                    break
                end
            end
            # println("   μ=$μ, total = $total")
            A = 2*pi*im*log(complex(z))^(s-1)/gamma(n)
            if  isreal(z) && real(z)>=1
                total -= A
            end
            # println("   μ=$μ, total = $total")
        elseif n==1
            total = -log(complex(1-z))
        elseif n==0
            total = z / (1-z)
        elseif n==-1
            total = z / (1-z)^2
        elseif n<-1
            # Crandall's 1.5 for s integer
            total = factorial(-n) * (-μ)^(n-1)
            tmp = 1
            for k=0:L
                # total -= μ^k * bernoulli(k-n+1, 0.0) / ( gamma(k+1)*(k-n+1) )
                total -= tmp * bernoulli(k-n+1, 0.0) / (k-n+1)
                tmp *= μ/(k+1)
                if abs(tmp)/abs(total) < 1.0e-30
                    break
                end
           end
        else
            error("Should not get this case")
        end

    # should have a case in here for when s is close to a real, positive integer
    # seed Wood 9.4
    # elseif abs(s - round(real(s))) < 0.01 # not sure where the right cut off is


    # could also use Wood 14.1 (square formula) to extend the range over which we sum directly - maybe???

    else
        # equivalent of Crandalls 1.4 for s non-integer, Wood (9.3)
        total = gamma(1-s) * (-μ)^(s-1)
        # println("   μ=$μ, total = $total")
        tmp = 1
        for k=0:L
            # total += μ^k * zeta(s-k)/factorial(Float64(k))
            total += tmp * zeta(s-k)
            # println("      tmp=$(tmp* zeta(s-k)),  total = $total")
            tmp *= μ/(k+1)
            if abs(tmp)/abs(total) < 1.0e-30
                break
            end
        end

        A = 2*pi*im*log(complex(z))^(s-1)/(gamma(s))
        if isreal(z) && real(z)>=1
            total -= A
        end
    end
    return total
end
