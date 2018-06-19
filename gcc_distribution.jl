struct DiscreteProba{T <: Real}
    func::Function
    values::Vector{T}
end

DiscreteProba(func::Function) = DiscreteProba(func, Float64[])

function (p::DiscreteProba)(k::Int)
    L = length(p.values)
    if k > L
        for j in max(1, L+1):k
            push!(p.values, p.func(j))
        end
    end
    return p.values[k]
end

function poisson(c, k)
    p = 1
    for j in 1:k
        p *= c/j
    end
    return exp(-c)*p
end

function gccpoisson(c, k)
    u = 0
    for _ in 1:1000
        u = exp(c*(u - 1))
    end
    return poisson(c, k)*(1 - u^k)
end

# Coeffs must be given from last to first
struct G1{T <: Real}
    coeffs::Vector{T}
end

function (P::G1)(z)
    poly = 0
    for (k, pk) in enumerate(P.coeffs)
        poly += k*pk*z^(k-1)
    end
    return poly
end

function ugen(r, u, tol=1e-15)
    Z = 0
    pk = 1
    k = 1
    uk_1 = 1  # u^(k-1)
    coeffs = Float64[]
    while k*pk > tol
        pk = r(k)/(1 - uk_1*u)
        push!(coeffs, pk)
        Z += k*pk
        uk_1 *= u
        k += 1
    end
    return G1(coeffs/Z)
end

function find_global_dist(r, u, tol=1e-10)
    ustart = -1
    res = nothing

    while abs(u - ustart) > tol
        ustart = u
        uold = -1
        ufunc = ugen(r, u)

        while abs(u - uold) > tol
            println(u)
            uold = u
            u = ufunc(u)
        end

        res = ufunc.coeffs
    end

    return u, res
end
