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

function single_param_SF(Λ0=Interval(nextfloat(2.), 2.5), Ls=[2, 3, 4], tol=0.005)
    for L in Ls
        generate_single_param_data(ScaleFreeGraph, Λ0, L, tol)
    end
end

function parse_root_data(data)
    Z, Λ = convert.(Vector{Float64}, data["bounds"])
    return Interval(Z...), Interval(Λ...), Symbol(data["status"])
end

function convert_data_to_S(graphtype, Ls)
    for L in Ls
        sp = SingleParameterGraph{graphtype}(L)
        Sdata = []
        data = JSON.parsefile("Plot generation/single_param_multiplex/$graphtype$L.json")

        for d in data
            Z, Λ, status = parse_root_data(d)
            push!(Sdata, Root(S(sp, Z, Λ) × Λ, status))
        end
        open("Plot generation/single_param_multiplex/$graphtype$(L)_S.json", "w") do file
            JSON.print(file, Sdata)
        end
    end
end
