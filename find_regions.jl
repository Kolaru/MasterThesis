if !isdefined(:first_run_)
    first_run_ = false
    include("main.jl")
end

using StaticArrays
using Plots
using NLsolve

const UNIT = IntervalBox(0..1, 2)

cut_corner(X) = bisect(bisect(X, 0.1)[1], 0.1)[1]

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

    f1(u1, u2, c1, c2) = 1 - (1 - g11(u1, c1))*(1 - g02(u2, c2))
    f2(u1, u2, c1, c2) = 1 - (1 - g01(u1, c1))*(1 - g12(u2, c2))
    F(U, C) = SVector(f1(U..., C...), f2(U..., C...))

    function detJ(u1, u2, c1, c2)
        j11 = 1 - dg11(u1, c1)*(1 - g02(u2, c2))
        j12 = - (1 - g11(u1, c1))*dg02(u2, c2)
        j21 = - dg01(u1, c1)*(1 - g12(u2, c2))
        j11 = 1 - (1 - g01(u1, c1))*dg12(u2, c2)
        return det([j11 j12 ; j21 j11])
    end

    function B!(res, u1, u2, c1, c2)
        res[1:2] = [u1, u2] - F([u1, u2], [c1, c2])
        res[3] = detJ(u1, u2, c1, c2)
    end

    numres = zeros(length(x_parameters))
    for (i, c) in enumerate(x_parameters)
        sol = nlsolve((res, UC) -> B!(res, UC..., c), [0, 0, 2.])
        numres[i] = sol.zero[3]
    end
    return numres
end

function find_region_in_two_layers(dist1, dist2, initial_region)
    ϵ = 1e-10  # Absolute tolerance for root search
    tol = 1e-2/2  # Minimal region diameter

    g01 = g0(dist1)
    g11 = g1(dist1)
    dg01 = dg0(dist1)
    dg11 = dg1(dist1)

    g02 = g0(dist2)
    g12 = g1(dist2)
    dg02 = dg0(dist2)
    dg12 = dg1(dist2)

    f1(u1, u2, c1, c2) = 1 - (1 - g11(u1, c1))*(1 - g02(u2, c2))
    f2(u1, u2, c1, c2) = 1 - (1 - g01(u1, c1))*(1 - g12(u2, c2))
    F(U, C) = SVector(f1(U..., C...), f2(U..., C...))

    T = typeof(initial_region)
    working = [initial_region]
    empties = T[]
    unkown =T[]
    sols = T[]

    while !isempty(working)
        C = shift!(working)

        if diam(C) < tol
            push!(unkown, C)
        else
            U = copy(UNIT)

            for i in 1:100
                U = F(U, C)
            end

            Utest = cut_corner(IntervalBox(U))
            if !any([X.lo > 1 - ϵ for X in Utest]) && IntervalBox(F(Utest, C)) ⊆ Utest
                push!(sols, C)
            elseif all([X.lo < 1 - ϵ for X in U])
                append!(working, bisect(C))
            else
                push!(empties, C)
            end
        end

        if length(working) > 10000
            warn("Limit on number of interval reached")
            append!(unkown, working)
            break
        end
    end
    return empties, working, sols
end
