using Espresso
import IterTools: product
import Base: widen

 function widen(X::Interval, relerr)
     lo_mod = X.lo > 0 ? 1 - relerr : 1 + relerr
     hi_mod = X.hi > 0 ? 1 + relerr : 1 - relerr

     return widen(Interval(lo_mod*X.lo, hi_mod*X.hi))
 end
widen(X::IntervalBox) = widen.(X)
widen(X::IntervalBox, relerr) = widen.(X, relerr)

abstract type Monotonicity end
struct Increasing <: Monotonicity end
struct Decreasing <: Monotonicity end

low_arg(mon::Type{Increasing}, X::Interval) = X.lo
low_arg(mon::Type{Decreasing}, X::Interval) = X.hi

high_arg(mon::Type{Increasing}, X::Interval) = X.hi
high_arg(mon::Type{Decreasing}, X::Interval) = X.lo

"""
    raw_extension(func::Function, Xs, monots)

Perform computation with the given argumentwise monotonic function `func` extended
for intervals. Arguments are specified as `(X, dir)`, where `X` is an input
interval and `dir` a type representing the direction of monoticity, either
`Increasing` or `Decreasing`.
"""
function raw_extension(func::Function, Xs, monots)
    low_args = (low_arg(mon, X) for (mon, X) in zip(monots, Xs))
    high_args = (high_arg(mon, X) for (mon, X) in zip(monots, Xs))

    low_bound = func(low_args...)
    high_bound = func(high_args...)

    # Conversion to Real is performed as some special function always return Complex
    try
        low_bound = convert(Real, low_bound)
        high_bound = convert(Real, high_bound)
    catch err
        if isa(err, InexactError)
            warn("Low bound: $low_bound")
            warn("High bound: $high_bound")
        end
        rethrow(err)
    end

    return Interval{Float64}(low_bound, high_bound)
end

function guaranteed_extension(func::Function, relerr, clampto,
            args::Vararg{Tuple{Interval{Float64}, Type{M} where M <: Monotonicity}, N}) where {N}

    Xs = (arg[1] for arg in args)
    monots = (arg[2] for arg in args)
    raw_res = raw_extension(func, Xs, monots)
    raw_res = widen(raw_res, relerr)

    clampto == -∞..∞ && return raw_res

    n = length(raw_res)
    clamped_intervals = Vector{Interval}(n)

    for i in 1:length(raw_res)
        lo_bound = max(raw_res[i].lo, clampto[i].lo)
        hi_bound = min(raw_res[i].hi, clampto[i].hi)
        clamped_intervals[i] = Interval(min(clampto[i].hi, lo_bound), max(clampto[i].lo, hi_bound))
    end

    if n == 1
        return clamped_intervals[1]
    else
        return IntervalBox(clamped_intervals...)
    end
end


function generate_extended_functions_expr(func, monots, relerr, clampto)
    argtype = Interval{Float64}
    N = length(monots)
    gen_args = [Symbol("X$i") for i in 1:N]

    typed_args = [subs(:(A::T), A=arg, T=argtype) for arg in gen_args]
    info_args = [subs(:((A, M)), A=arg, M=monot) for (arg, monot) in zip(gen_args, monots)]

    func_def = :( FUNC(ARGS...) = guaranteed_extension(FUNC, RELERR, CLAMPTO, INFOS...) )
    func_expr = subs(func_def, FUNC=func, ARGS=typed_args, INFOS=info_args,
                               RELERR=relerr, CLAMPTO=clampto)

    intervalled_args = [subs(:(Interval(A)), A=arg) for arg in gen_args]
    both_types = zip(typed_args, gen_args)
    type_combinations = collect(product(both_types...))

    partial_def = :( FUNC(ARGS...) = FUNC(I_ARGS...) )
    partial_expr = Expr[]

    for combination in type_combinations[2:(end-1)]
        combination = collect(combination)
        pexpr = subs(partial_def, FUNC=func, ARGS=combination, I_ARGS=intervalled_args)
        push!(partial_expr, pexpr)
    end

    return quote
        $func_expr
        $(partial_expr...)
    end
end


function infer_monotonicity(f, domain)
    N = length(domain)

    mids = mid(domain)
    los = [dom.lo for dom in domain]
    his = [dom.hi for dom in domain]
    monots = Array{Type}(N)

    for i in 1:N
        lo_args = collect(mids)
        hi_args = collect(mids)
        lo_args[i] = los[i]
        hi_args[i] = his[i]

        if f(lo_args...) < f(hi_args...)
            monots[i] = Increasing
        else
            monots[i] = Decreasing
        end
    end

    return monots
end


"""
    @extend_monotonic f(dir [,...]) clampto=C domain=D

Define `2^N-1` new funcitons (where `N` is the number of given arguments of `f`)
extending the argumentwise monotonic function `f` for intervals computations for
any combination of `Real` and `Interval{Float64}` arguments.

Arguments `dir` of `f` are the direction of the monoticity for that argument,
either `Increasing` or `Decreasing`.
"""
macro extend_monotonic(func, domain, relerr=:(0), clampto=:(-∞..∞))
    func = esc(func)
    domain = eval(domain)
    N = length(domain)
    relerr = eval(relerr)
    clampto = eval(clampto)

    monots = infer_monotonicity(f, domain)

    return generate_extended_functions_expr(func, monots, relerr, clampto)
end
