using StaticArrays
using IntervalArithmetic
using IntervalRootFinding
using IterTools

import IntervalRootFinding: BreadthFirstSearch
import JSON
import ForwardDiff: jacobian

using GeneratingFunctions
using Graphs

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
    return SVector{L}(1. .- res)
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

function ufunc_single_param(z, layer, param, L)
    return 1. - (1. - g1(layer, z, param))*(1 - g0(layer, z, param))^(L - 1)
end

function ures_single_param(z, layer, param, L)
    return ufunc_single_param(z, layer, param, L) - z
end

function refine_single_param(R, layer, L, level=6)
    working = [R]
    final = []
    k = 1
    while k < 2^level && !isempty(working)
        k += 1
        X = interval(popfirst!(working))
        U, P = X
        CU = ufunc_single_param(U, layer, P, L)

        if CU ⊆ U
            push!(final, Root(X, :exist))
        elseif CU ∩ U != ∅
            X1, X2 = bisect(Root(X, :unknown))
            push!(working, X1)
            push!(working, X2)
        end
    end

    append!(final, working)
    return final
end

function JSON.lower(rt::Root)
    X = interval(rt)
    N = length(X)
    bounds = zeros(2, N)
    for (k, x) in enumerate(X)
        bounds[1, k] = x.lo
        bounds[2, k] = x.hi
    end
    return Dict("bounds" => bounds, "status" => rt.status)
end

function generate_single_param_data(layer, pbounds, L=2, tol=0.01)
    function f(zp)
        res = ures_single_param(zp[1], layer, zp[2], L)
        return SVector{2}(fill(res, 2))
    end
    R = UNITINTERVAL × pbounds
    C = Krawczyk(f, zp -> jacobian(f, zp))

    # A search takes a starting region, a contractor and a tolerance as argument
    search = BreadthFirstSearch(R, C, tol)

    global endtree
    for (k, tree) in enumerate(search)
        # println(tree)  # The tree has custom printing
        global endtree = tree
    end
    rts = data(endtree)
    refined_rts = []
    for (k, rt) in enumerate(rts)
        refined = refine_single_param(rt, layer, L)
        append!(refined_rts, refined)
    end

    open("Plot generation/single_param_multiplex/$layer$L.json", "w") do f
        JSON.print(f, refined_rts)
    end
    return rts

end
