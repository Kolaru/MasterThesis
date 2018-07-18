using Espresso
import IterTools: product
drop = Iterators.drop

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

    N = length(expr_dict[:monots])
    func = expr_dict[:func]

    gen_args = [Symbol("X$i") for i in 1:N]
    argtype = Interval{Float64}

    typed_args = [subs(:(A::T), A=arg, T=argtype) for arg in gen_args]
    info_args = [subs(:((A, M)), A=arg, M=monot) for (arg, monot) in zip(gen_args, expr_dict[:monots])]

    Ifunc_expr = quote
        function $(esc(func))($(typed_args...))
            monotonic_extender($(esc(func)), $(info_args...))
        end
    end

    typed_reals = [subs(:(A::Real), A=arg) for arg in gen_args]
    intervalled_args = [subs(:(Interval(A)), A=arg) for arg in gen_args]
    both_types = [[typed_args[i], typed_reals[i]] for i in 1:N]

    type_combinations = collect(product(both_types...))

    partial_expr = Expr[]

    for combination in type_combinations[2:(end-1)]
        pexpr = quote
            function $(esc(func))($(combination...))
                 $(esc(func))($(intervalled_args...))
            end

        end
        push!(partial_expr, pexpr)
    end

    return quote
        $Ifunc_expr

        $(partial_expr...)
    end
end
