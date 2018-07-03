if !isdefined(:first_run_)
    first_run_ = false
    include("main.jl")
end

using StaticArrays
using Plots

rect(X::IntervalBox) = rect(X...)
rect(W, H) = Shape([W.lo, W.lo, W.hi, W.hi], [H.lo, H.hi, H.hi, H.lo])
rectangle(w, h, x, y) = Shape(x + [0,w,w,0], y + [0,0,h,h])

ϵ = 1e-5
tol = 1e-2

g = g0(ErdosRenyiGraph)

f1(u1, u2, c1, c2) = 1 - (1 - g(u1, c1))*(1 - g(u2, c2))
f2(u1, u2, c1, c2) = 1 - (1 - g(u1, c1))*(1 - g(u2, c2))
F(U, C) = SVector(f1(U..., C...), f2(U..., C...))

UNIT = IntervalBox(0..1, 2)

Cinit = IntervalBox(1..3, 2)
working = [Cinit]
results = []

while !isempty(working)
    C = shift!(working)

    if diam(C) < tol
        push!(results, C)
    else
        U = copy(UNIT)

        for i in 1:100
            U = F(U, C)
        end

        if all([X.lo < 1 - ϵ for X in U])
            append!(working, bisect(C))
        end
    end

    if length(working) > 10000
        break
    end
end
