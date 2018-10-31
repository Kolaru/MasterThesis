using IntervalArithmetic
using IntervalRootFinding
using IterTools

import JSON

using GeneratingFunctions
using Graphs
using Utils

"""
    Fixpoint function `ψ`.

Order of arguments `λ` and `z` are inverted compare to what is done in the
report, to match the definition of `g0` and `g1` in `GeneratingFunctions.jl`.
"""
function ψ(graphtype, z, λ, L)
    return 1. - (1. - g1(graphtype, z, λ))*(1 - g0(graphtype, z, λ))^(L - 1)
end

function dψ(graphtype, z, λ, L)
    p0 = 1 - g0(graphtype, z, λ)
    return ( dg1(graphtype, z, λ)*p0 + (L - 1)*(1 - g1(graphtype, z, λ))*dg0(graphtype, z, λ) ) * p0^(L-2)
end

function fixpoint_refine(graphtype, X, Λ, L)
    FX = ψ(graphtype, X, Λ, L)
    status = :unkown
    if FX ⊆ X
        satus = :exist
    end

    NX = FX ∩ X
    return Root(NX × Λ, status)
end

"""
    Residual function `f_λ`.

`f_λ(z) == 0` iff `ψ(z, λ) = z`.
"""
struct ResFunc{GT <: GraphType}
    graphtype::Type{GT}
    λ
    L::Int
end

func(rf::ResFunc) = z -> ψ(rf.graphtype, z, rf.λ, rf.L) - z
deriv(rf::ResFunc) = z -> dψ(rf.graphtype, z, rf.λ, rf.L) - 1


function generate_single_param_data(graphtype, Λ0, L=2, tol=0.005)
    working = [UNITINTERVAL × Λ0]
    stored = []

    while !isempty(working)
        X, Λ = popfirst!(working)
        resfunc = ResFunc(graphtype, Λ, L)
        f = func(resfunc)
        df = deriv(resfunc)

        rts = roots(f, df, X, Krawczyk, tol)

        for rt in rts
            Y = rt.interval
            if isunique(rt)
                push!(stored, Root(Y × Λ, :unique))
            else
                if diam(Y × Λ) < tol
                    refined = fixpoint_refine(graphtype, Y, Λ, L)
                    !isempty(refined.interval) && push!(stored, refined)
                else
                    append!(working, bisect(Y × Λ))
                end
            end
        end
    end

    open("Plot generation/single_param_multiplex/$graphtype$L.json", "w") do file
        JSON.print(file, stored)
    end

    return stored
end

function single_param_ER(Λ0=2..4, Ls=[2, 3, 4, 5], tol=0.02)
    for L in Ls
        generate_single_param_data(ErdosRenyiGraph, Λ0, L, tol)
    end
end

function single_param_geometric(Λ0=2..5, Ls=[2, 3, 4], tol=0.02)
    for L in Ls
        generate_single_param_data(GeometricGraph, Λ0, L, tol)
    end
end

function single_param_SF(Λ0=1.5..2.4, Ls=[2, 3, 4], tol=0.005)
    for L in Ls
        generate_single_param_data(ScaleFreeGraph, Λ0, L, tol)
    end
end
