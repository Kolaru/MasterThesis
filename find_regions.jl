using IntervalArithmetic
using IntervalRootFinding
using NLsolve
using StaticArrays

import ForwardDiff: jacobian
import JSON
import JSON: json

using GeneratingFunctions
using Graphs
using Utils

rect(X::IntervalBox) = rect(X...)
rect(W, H) = Shape([W.lo, W.lo, W.hi, W.hi], [H.lo, H.hi, H.hi, H.lo])

function numerical_solve_critical_curve(dist1, dist2, x_parameters)
    g01 = g0(dist1)
    g11 = g1(dist1)
    dg01 = dg0(dist1)
    dg11 = dg1(dist1)

    g02 = g0(dist2)
    g12 = g1(dist2)
    dg02 = dg0(dist2)
    dg12 = dg1(dist2)

    z1func(z1, z2, λ1, λ2) = 1 - (1 - g11(z1, λ1))*(1 - g02(z2, λ2))
    z2func(z1, z2, λ1, λ2) = 1 - (1 - g01(z1, λ1))*(1 - g12(z2, λ2))
    Ufunc(U, C) = SVector(z1func(U..., C...), z2func(U..., C...))

    function detJ(z1, z2, λ1, λ2)
        j11 = 1 - dg11(z1, λ1)*(1 - g02(z2, λ2))
        j12 = - (1 - g11(z1, λ1))*dg02(z2, λ2)
        j21 = - dg01(z1, λ1)*(1 - g12(z2, λ2))
        j11 = 1 - (1 - g01(z1, λ1))*dg12(z2, λ2)
        return det([j11 j12 ; j21 j11])
    end

    function B!(res, z1, z2, λ1, λ2)
        res[1:2] = [z1, z2] - Ufunc([z1, z2], [λ1, λ2])
        res[3] = detJ(z1, z2, λ1, λ2)
    end

    numres = zeros(length(x_parameters))
    for (i, c) in enumerate(x_parameters)
        sol = nlsolve((res, UC) -> B!(res, UC..., c), [0, 0, 2.])
        numres[i] = sol.zero[3]
    end
    return numres
end

contains_trivial(X::IntervalBox{N, T}) where {N, T} = all(ones(N) .∈ X)

function nontrivial_sol(X::IntervalBox{N, T}, status) where {N, T}
    return status == :exist && !contains_trivial(X)
end

approx_trivial(X, ε) = all([x.lo > 1 - ε for x in X])
is_trivial(rt, ε) = isunique(rt) && contains_trivial(rt)

function Base.isapprox(X::Interval{T}, Y::Interval{T} ;
                  rtol::Real=√eps(T), atol::Real=zero(T)) where {T <: Real}
    #
    maxerr = max(atol, rtol*max(mag(X), mag(Y)))
    return abs(X.lo - Y.lo) + abs(X.hi - Y.hi) <= maxerr
end

function Base.isapprox(X::IntervalBox{N, T}, Y::IntervalBox{N, T} ;
                  rtol::Real=√eps(T), atol::Real=zero(T)) where {N, T <: Real}
    #
    return all(X .≈ Y)
end

function psi_process(ψ, X0, Λ, existstop=true)
    prevX = X0
    X = ψ(X0, Λ)
    status = :unkown

    while !(prevX ≈ X)
        prevX = X
        X = ψ(X, Λ)

        if X ⊆ prevX
            status = :exist
            existstop && break
        end
    end

    return X, status
end

function find_region(ψ, Λ0::IntervalBox{N, TT}, δ, ε) where {N, TT}
    working = [Λ0]
    T = IntervalBox{N, TT}[]
    NT = IntervalBox{N, TT}[]
    U = IntervalBox{N, TT}[]

    while !isempty(working)
        print("\rT : $(length(T))   NT : $(length(NT))   U : $(length(U))   W : $(length(working))  ")
        Λ = pop!(working)

        if diam(Λ) < δ
            push!(U, Λ)
            continue
        end

        X, status = psi_process(ψ, IntervalBox(UNITINTERVAL, N), Λ, false)
        if approx_trivial(X, ε)
            push!(T, Λ)
            continue
        end

        # We use the lower corner of the box
        Y = IntervalBox([x.lo..(0.25 + 0.75*x.lo) for x in X])
        CY, status = psi_process(ψ, Y, Λ)

        if nontrivial_sol(CY, status)
            push!(NT, Λ)
        else
            append!(working, bisect(Λ))
        end
    end

    println()
    return T, NT, U
end

function find_region_in_two_layers(dist1, dist2, Λ0, δ=5e-3, ε=1e-3)
    g01 = g0(dist1)
    g11 = g1(dist1)

    g02 = g0(dist2)
    g12 = g1(dist2)

    T, NT, U = find_region(ψ([dist1, dist2]), Λ0, δ, ε)

    println("T : $(length(T))")
    println("NT : $(length(NT))")
    println("U : $(length(U))")

    data = Dict("unkown" => U,
                "trivial" => T,
                "nontrivial" => NT)

    open("Plot generation/critical_region/$(dist1)_$dist2.json", "w") do file
        JSON.print(file, data)
    end

    return U, T, NT
end

function find_region_in_monolayer(dist, initial_region)
    Ufunc(u, C) = g1(dist, u, C[1], C[2])

    return find_region(Ufunc, initial_region, Val{1})
end
