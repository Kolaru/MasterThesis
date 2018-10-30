using IntervalArithmetic
using IntervalRootFinding
using IterTools

import JSON

using GeneratingFunctions
using Graphs
using Utils

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
