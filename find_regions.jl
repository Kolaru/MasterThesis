if !isdefined(:first_run_)
    first_run_ = false
    include("main.jl")
end

using StaticArrays
using Plots
using NLsolve

rect(X::IntervalBox) = rect(X...)
rect(W, H) = Shape([W.lo, W.lo, W.hi, W.hi], [H.lo, H.hi, H.hi, H.lo])
rectangle(w, h, x, y) = Shape(x + [0,w,w,0], y + [0,0,h,h])

ϵ = 1e-10
tol = 1e-2/2

g = g0(ErdosRenyiGraph)
dg = dg0(ErdosRenyiGraph)

f1(u1, u2, c1, c2) = 1 - (1 - g(u1, c1))*(1 - g(u2, c2))
f2(u1, u2, c1, c2) = 1 - (1 - g(u1, c1))*(1 - g(u2, c2))
F(U, C) = SVector(f1(U..., C...), f2(U..., C...))

function detJ(u1, u2, c1, c2)
    j11 = 1 - dg(u1, c1)*(1 - g(u2, c2))
    j12 = - dg(u2, c2)*(1 - g(u1, c1))
    j21 = - dg(u1, c1)*(1 - g(u2, c2))
    j11 = 1 - dg(u2, c2)*(1 - g(u1, c1))
    return det([j11 j12 ; j21 j11])
end

function B!(res, u1, u2, c1, c2)
    res[1:2] = [u1, u2] - F([u1, u2], [c1, c2])
    res[3] = detJ(u1, u2, c1, c2)
end

cc = linspace(2.1, 3, 100)
numres = zeros(100)
for (i, c) in enumerate(cc)
    sol = nlsolve((res, UC) -> B!(res, UC..., c), [0, 0, 2.], autodiff=true)
    numres[i] = sol.zero[3]
end

UNIT = IntervalBox(0..1, 2)

Cinit = IntervalBox(1..3, 2)
working = [Cinit]
unkown = []
sols = []

halve(X) = bisect(bisect(X, 0.1)[1], 0.1)[1]
isin(A, B) = all([b.lo <= a.lo && a.hi <= b.hi for (a, b) in zip(A, B)])

while !isempty(working)
    C = shift!(working)

    if diam(C) < tol
        push!(unkown, C)
    else
        U = copy(UNIT)

        for i in 1:100
            U = F(U, C)
        end

        Utest = halve(IntervalBox(U))
        if !any([X.lo > 1 - ϵ for X in Utest]) && IntervalBox(F(Utest, C)) ⊆ Utest
            push!(sols, C)
        elseif all([X.lo < 1 - ϵ for X in U])
            append!(working, bisect(C))
        end
    end

    if length(working) > 10000
        break
    end
end
