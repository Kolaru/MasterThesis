using Espresso
import IterTools: product

abstract type Monotonicity end
struct Increasing <: Monotonicity end
struct Decreasing <: Monotonicity end

low_arg(mon::Type{Increasing}, X::Interval) = X.lo
low_arg(mon::Type{Decreasing}, X::Interval) = X.hi

high_arg(mon::Type{Increasing}, X::Interval) = X.hi
high_arg(mon::Type{Decreasing}, X::Interval) = X.lo

"""
    monotonic_extender(func, arginfos::Vararg{Tuple{Interval{Float64}, Symbol}, N}) where {N}

Perform computation with the given argumentwise monotonic function `func` extended
for intervals. Arguments are specified as `(X, dir)`, where `X` is an input
interval and `dir` a type representing the direction of monoticity, either
`Increasing` or `Decreasing`.
"""
function monotonic_extender(func::Function, args::Vararg{Tuple{Interval{Float64}, Type{M} where M <: Monotonicity}, N}) where {N}
    Xs = (Interval{BigFloat}(arg[1]) for arg in args)
    monotonicities = (arg[2] for arg in args)

    low_args = (low_arg(mon, X) for (mon, X) in zip(monotonicities, Xs))
    high_args = (high_arg(mon, X) for (mon, X) in zip(monotonicities, Xs))

    # Conversion to Real is performed as some special function always return Complex
    low_bound = convert(Real, func(low_args...))
    high_bound = convert(Real, func(high_args...))

    return Interval{Float64}(low_bound, high_bound)
end

function generate_extended_functions_expr(func, monots)
    argtype = Interval{Float64}
    N = length(monots)
    gen_args = [Symbol("X$i") for i in 1:N]

    typed_args = [subs(:(A::T), A=arg, T=argtype) for arg in gen_args]
    info_args = [subs(:((A, M)), A=arg, M=monot) for (arg, monot) in zip(gen_args, monots)]

    func_def = :( FUNC(ARGS...) = monotonic_extender(FUNC, INFOS...) )
    func_expr = subs(func_def, FUNC=func, ARGS=typed_args, INFOS=info_args)

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

"""
    @extend_monotonic f(dir [,...])

Define `2^N-1` new funcitons (where `N` is the number of given arguments of `f`)
extending the argumentwise monotonic function `f` for intervals computations for
any combination of `Real` and `Interval{Float64}` arguments.

Arguments `dir` of `f` are the direction of the monoticity for that argument,
either `Increasing` or `Decreasing`.
"""
macro extend_monotonic(expr)
    expr_dict = matchex(:(func(monots...)), expr, phs=[:func, :monots])
    expr_dict = get(expr_dict)

    func = esc(expr_dict[:func])
    monots = expr_dict[:monots]

    return generate_extended_functions_expr(func, monots)
end
