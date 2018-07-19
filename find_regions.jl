using StaticArrays
using Plots
using NLsolve

const UNIT = 0..1

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

    u1func(u1, u2, c1, c2) = 1 - (1 - g11(u1, c1))*(1 - g02(u2, c2))
    u2func(u1, u2, c1, c2) = 1 - (1 - g01(u1, c1))*(1 - g12(u2, c2))
    Ufunc(U, C) = SVector(u1func(U..., C...), u2func(U..., C...))

    function detJ(u1, u2, c1, c2)
        j11 = 1 - dg11(u1, c1)*(1 - g02(u2, c2))
        j12 = - (1 - g11(u1, c1))*dg02(u2, c2)
        j21 = - dg01(u1, c1)*(1 - g12(u2, c2))
        j11 = 1 - (1 - g01(u1, c1))*dg12(u2, c2)
        return det([j11 j12 ; j21 j11])
    end

    function B!(res, u1, u2, c1, c2)
        res[1:2] = [u1, u2] - Ufunc([u1, u2], [c1, c2])
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
    g01 = g0(dist1)
    g11 = g1(dist1)

    g02 = g0(dist2)
    g12 = g1(dist2)

    u1func(u1, u2, c1, c2) = 1 - (1 - g11(u1, c1))*(1 - g02(u2, c2))
    u2func(u1, u2, c1, c2) = 1 - (1 - g01(u1, c1))*(1 - g12(u2, c2))
    Ufunc(U, C) = SVector(u1func(U..., C...), u2func(U..., C...))

    return find_region(Ufunc, initial_region)
end

function find_region_in_monolayer(dist, initial_region)
    Ufunc(u, C) = g1(dist, u, C[1], C[2])

    return find_region(Ufunc, initial_region, 1)
end

function find_region(Ufunc, initial_region, Udim=2)
    ϵ = 1e-10  # Absolute tolerance for root search
    tol = 1e-2/2  # Minimal region diameter

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
            if Udim == 1
                U = deepcopy(UNIT)
            else
                U = IntervalBox(UNIT, Udim)
            end

            for i in 1:100
                U = Ufunc(U, C)
            end

            Utest = cut_corner(IntervalBox(U))
            if !any([X.lo > 1 - ϵ for X in Utest]) && IntervalBox(Ufunc(Utest, C)) ⊆ Utest
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
