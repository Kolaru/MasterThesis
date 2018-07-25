module PowerLawDistribution

import SpecialFunctions: zeta

export plrand

"""
    Compute CDF for a pure power-law distributed discrete random variable.

    Take a list of already calculated values of the zeta function as argument
    to avoid recomputing them.
"""
function Ppl(x, α, known_zetas=Dict{Int, Float64}())
    k = keys(known_zetas)
    if !(1 in k)
        known_zetas[1] = zeta(α, 1)
    end

    if !(x in k)
        known_zetas[x] = zeta(α, x)
    end

    return known_zetas[x]/known_zetas[1]
end


"""
    Generate `n` random power-law distributed integers.

    `α` must be bigger than 1.
    Algorithm from Clauset 2009 (Appendix D).
"""
function plrand(α, n)
    if α <= 1
        error("Power-law generator : α smaller than 1 causes non converging normalization constant.")
    end

    rs = rand(n)
    res = zeros(Int, n)
    known_zetas = Dict{Int, Float64}()

    for (i, r) = enumerate(rs)
        low = 1
        high = 2*low
        while Ppl(high, α, known_zetas) >= r
            low = high
            high *= 2
        end

        while high - low > 1 && low > 0
            mid = low + div(high-low, 2)
            if Ppl(mid, α, known_zetas) < r
                high = mid
            else
                low = mid
            end
        end
        res[i] = low
    end

    return res
end

end
