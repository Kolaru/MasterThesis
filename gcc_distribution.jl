import StatsBase: fit, Histogram

struct DiscreteProba{F <: Function, T <: Real}
    func::F
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
    p = 1.
    for j in 1:k
        p *= c/j
    end
    return exp(-c)*p
end

function gccpoisson(c, k)
    u = 0.
    for _ in 1:1000
        u = exp(c*(u - 1))
    end
    return poisson(c, k)*(1 - u^k)
end

# Coeffs must be given from first to last
struct G1{T <: Real}
    coeffs::Vector{T}
end

function (P::G1)(z)
    poly = 0.
    for (k, pk) in reverse(collect(enumerate(P.coeffs)))
        poly += k*pk*z^(k-1)
    end
    return poly
end

function ugen(r::R, u, N) where {R <: Function}
    if u == 0.
        return G1([r(k) for k in 1:N])
    end
    Z = 0.
    pk = 1.
    uk_1 = u^(N-1)  # u^(k-1)
    coeffs = Float64[]
    for k in N:-1:1
        pk = r(k)/(1 - uk_1*u)
        push!(coeffs, pk)
        Z += k*pk
        uk_1 /= u
    end
    return G1(reverse(coeffs)/Z)
end

function find_global_dist(r::R, u=0.0, tol=1e-14, N=100) where {R <: Function}
    ustart = -1.
    res = Float64[]
    ufunc = ugen(r, u, N)

    while abs(u - ustart)/u > tol
        firstpass = false
        ustart = u
        u = ufunc(u)
        ufunc = ugen(r, u, N)
        res = ufunc.coeffs
    end

    return u, res/sum(res)
end

function find_global_dist(r::AbstractArray, u=0.0, tol=1e-14, N=100)
    return find_global_dist(k -> r[k], u, tol, N)
end

function mimic_real_network(g::Graph)
    deg = degrees(g)
    hist = fit(Histogram, deg, nbins=maximum(deg), closed=:left)
    r = normalize(hist).weights
    return find_global_dist(r, 0.0, 1e-14, length(r))
end
