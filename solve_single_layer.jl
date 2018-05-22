push!(LOAD_PATH,pwd())

using ValidatedNumerics
import IntervalRootFinding.clean_roots

using GeneratingFunctions

const UNITINTERVAL = 0..1
const TOL = 1e-12

low(interval) = interval.lo
high(interval) = interval.hi


function find_u(dist)
    X = UNITINTERVAL
    rts = roots(u -> u - dist.g1(u), X, Newton(u -> 1 - dist.dg1(u)), TOL)
    return rts
end

function find_S(dist)
    uu = find_u(dist)
    corr_u = []
    for u in uu
        if !is_unique(u)
            if 1.0 in u.interval
                push!(corr_u, Root(u.interval, :unique))
            else
                warn("Unkown intervals for degree distribution $dist")
                @show u.interval
                println(1.0 in u.interval)
                break
            end
        end
    end
    uu = [u for u in uu if is_unique(u)]
    append!(uu, corr_u)

    if isempty(uu)
        return @interval 0
    else
        SS = map(z -> 1 - dist.g0(z.interval), uu)
        return SS
    end
end

info("Initialization done.")

function solve_on_range(Dist, cc)
    H = []
    L = []

    for c in cc
        println(c)
        dist = Dist(c)
        SS = find_S(dist)
        push!(H, maximum(SS))
        push!(L, minimum(SS))
    end
    return L, H
end
