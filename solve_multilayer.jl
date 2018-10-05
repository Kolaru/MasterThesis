using StaticArrays
using IntervalArithmetic
using IntervalRootFinding
using IterTools

include("network_generation.jl")
include("src/generating_functions/generating_functions.jl")

const UNITINTERVAL = 0..1

low(interval) = interval.lo
high(interval) = interval.hi

function ufunc(zz::T, layers, params) where T
    L = length(layers)
    res = ones(typeof(zz[1]), L)

    for (j, layer, param) in zip(1:L, layers, params)
        for i in 1:L
            if i == j
                res[j] *= 1. - g1(layer, zz[i], param...)
            else
                res[j] *= 1. - g0(layer, zz[i], param...)
            end
        end
    end
    return SVector{L}(1. - res)
end

function ures(zz::SVector, layers, params)
    L = length(layers)
    return zz .- ufunc(zz, layers, params)
end

function ures(zz::Region, layers, params)
    L = length(layers)
    X = zz .- ufunc(zz, layers, params)
    return X.v
end

function find_u(layers, params)
    L = length(layers)
    X = IntervalBox(UNITINTERVAL, L)
    rts = roots(uu -> ures(uu, layers, params), X, Krawczyk, 1e-3)
    return rts
end

function find_S(layers, params)
    L = length(layers)
    uu = find_u(layers, params)
    res = ones(typeof(uu[1]), L)

    for rt in uu
        for (layer, u) in zip(layers, rt.interval)
            res[i] *= 1. - g0(layer, uu[i], param...)
        end
    end

    return SS
end


function ures_single_param(z, layer, param, L)
    return 1. - z - (1. - g1(layer, z, param))*(1 - g0(layer, z, param))^(L - 1)
end

function generate_single_param_data(layer, pbounds, L=2, tol=0.5)
    function f(zp)
        res = ures_single_param(zp[1], layer, zp[2], L)
        println(zp)
        return SVector{2}(res, res)
    end
    R = UNITINTERVAL × pbounds
    println(f(R))
    rts = roots(f, R, Krawczyk, tol)
    return rts
end
