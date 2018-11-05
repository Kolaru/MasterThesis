using IntervalArithmetic
using IntervalRootFinding
using IterTools

import JSON

using GeneratingFunctions
using Graphs
using Utils

function fixpoint_refine(spgraph, X, Λ)
    FX = ψ(spgraph, X, Λ)
    status = :unkown
    if FX ⊆ X
        satus = :exist
    end

    NX = FX ∩ X
    return Root(NX × Λ, status)
end

function compute_data(graphtype, Λ0, L, tol)
    working = [UNITINTERVAL × Λ0]
    stored = []
    spgraph = SingleParameterGraph{graphtype}(L)

    while !isempty(working)
        X, Λ = popfirst!(working)
        f(z) = ψ(spgraph, z, Λ) - z
        df(z) = dψ(spgraph, z, Λ) - 1

        rts = roots(f, df, X, Krawczyk, tol)

        for rt in rts
            Y = rt.interval
            if isunique(rt)
                push!(stored, Root(Y × Λ, :unique))
            else
                if diam(Y × Λ) < tol
                    refined = fixpoint_refine(spgraph, Y, Λ)
                    !isempty(refined.interval) && push!(stored, refined)
                else
                    append!(working, bisect(Y × Λ))
                end
            end
        end
    end

    return stored
end

function generate_single_param_data(graphtype, Λ0, L=2, tol=0.005)
    stored = compute_data(graphtype, Λ0, L, tol)

    open("Plot generation/single_param_multiplex/$graphtype$L.json", "w") do file
        JSON.print(file, stored)
    end

    return stored
end

function generate_single_param_data(graphtype, Λ0s::Vector, L=2, tol=0.005)
    stored = []
    for Λ0 in Λ0s
        append!(stored, compute_data(graphtype, Λ0, L, tol))
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

function single_param_SF(Λ0=1.2..2.4, Ls=[2, 3, 4], tol=0.005)
    Λ0s = []
    if Λ0.lo < 2 && Λ0.hi > 2
        push!(Λ0s, Interval(Λ0.lo, prevfloat(2.0)))
        push!(Λ0s, Interval(nextfloat(2.0), Λ0.hi))
    else
        push!(Λ0s, Λ0)
    end

    for L in Ls
        data = generate_single_param_data(ScaleFreeGraph, Λ0s, L, tol)

        L == 3 && println(data)
    end
end
