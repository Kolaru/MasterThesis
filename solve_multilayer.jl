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

function dufunc_single_param(z, layer, param, L)
    p0 = 1 - g0(layer, z, param)
    return ( dg1(layer, z, param)*p0 + (L - 1)*(1 - g1(layer, z, param))*dg0(layer, z, param) ) * p0^(L-2)
end

function ures_single_param(z, layer, param, L)
    return ufunc_single_param(z, layer, param, L) - z
end

function dures_single_param(z, layer, param, L)
    return dufunc_single_param(z, layer, param, L) - 1
end

struct PhiContractor <: Contractor{Function}
    layer
    Λ
    L
end

function (C::PhiContractor)(r, tol)
    layer = C.layer
    Λ = C.Λ
    L = C.L

    X = interval(r)
    former_status = root_status(r)
    PX = ufunc_single_param(X, layer, Λ, L)

    NX = PX ∩ X

    isempty(NX) && return Root(X, :empty)
    isinf(X) && return Root(X, :unkown)  # force bisection

    if PX ⊆ X  # isinterior; know there's a unique root inside
        return Root(NX, :exist)
    end

    return Root(NX, :unkown)
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

function fixpoint_refine(X, layer, Λ, L)
    FX = ufunc_single_param(X, layer, Λ, L)
    status = :unkown
    if FX ⊆ X
        satus = :exist
    end

    NX = FX ∩ X
    return Root(NX × Λ, status)
end

function generate_single_param_data(layer, pbound, L=2, tol=0.005, C=Krawczyk)
    working = [UNITINTERVAL × pbound]
    stored = []

    while !isempty(working)
        X, Λ = popfirst!(working)
        f(z) = ures_single_param(z, layer, Λ, L)
        df(z) = dures_single_param(z, layer, Λ, L)

        rts = roots(f, df, X, Krawczyk, tol)

        for rt in rts
            Y = rt.interval
            if isunique(rt)
                push!(stored, Root(Y × Λ, :unique))
            else
                if diam(Y × Λ) < tol
                    refined = fixpoint_refine(Y, layer, Λ, L)
                    !isempty(refined.interval) && push!(stored, refined)
                else
                    append!(working, bisect(Y × Λ))
                end
            end
        end
    end

    open("Plot generation/single_param_multiplex/$layer$L.json", "w") do file
        JSON.print(file, stored)
    end

    return stored
end

function single_param_ER(pbound=2..4, Ls=[2, 3, 4, 5], tol=0.02)
    for L in Ls
        generate_single_param_data(ErdosRenyiGraph, pbound, L, tol)
    end
end

function single_param_geometric(pbound=2..5, Ls=[2, 3, 4], tol=0.02)
    for L in Ls
        generate_single_param_data(GeometricGraph, pbound, L, tol)
    end
end

function single_param_SF(pbound=1.5..2.4, Ls=[2, 3, 4], tol=0.005)
    for L in Ls
        println("L = $L")
        regions = generate_single_param_data(ScaleFreeGraph, pbound, L, tol)
        println(length(regions))
        println()
    end
end

function generate_single_param_data_SF(pbound, L=2, tol=0.005)
    working = [UNITINTERVAL × pbound]
    stored = []
    layer = ScaleFreeGraph

    while !isempty(working)
        X, Λ = popfirst!(working)
        C = PhiContractor(layer, Λ, L)

        search = BreadthFirstSearch(X, C, tol)

        global endtree
        for tree in search
            global endtree = tree
        end
        rts = data(endtree)

        for rt in rts
            Y = rt.interval
            if isunique(rt)
                push!(stored, Root(Y × Λ, :unique))
            else
                if diam(Y × Λ) < tol
                    refined = fixpoint_refine(Y, layer, Λ, L)
                    !isempty(refined.interval) && push!(stored, refined)
                else
                    append!(working, bisect(Y × Λ))
                end
            end
        end
    end

    open("Plot generation/single_param_multiplex/$layer$L.json", "w") do file
        JSON.print(file, stored)
    end

    return stored
end
