using Espresso
import IterTools: product

drop = Iterators.drop

lo(X::Interval) = X.lo
hi(X::Interval) = X.hi

function low_bound_func(info::Symbol)
    if info == :up
        return lo
    elseif info == :down
        return hi
    end
end

function high_bound_func(info::Symbol)
    if info == :up
        return hi
    elseif info == :down
        return lo
    end
end

"""
    monotonic_extender(func, arginfos::Vararg{Tuple{Interval{Float64}, Symbol}, N}) where {N}

Perform computation with the given elementwise monotonic function `func` extended
for intervals. Arguments are specified as `(X, dir)`, where `X` is an input interval and `dir` a symbol representing the direction of monoticity, either
`:up` or `:down`.
"""
function monotonic_extender(func, arginfos::Vararg{Tuple{Interval{Float64}, Symbol}, N}) where {N}

    low_funcs = (low_bound_func(arg[2]) for arg in arginfos)
    high_funcs = (high_bound_func(arg[2]) for arg in arginfos)

    Bargs = (Interval{BigFloat}(arg[1]) for arg in arginfos)
    low_args = (f(arg) for (f, arg) in zip(low_funcs, Bargs))
    high_args = (f(arg) for (f, arg) in zip(high_funcs, Bargs))

    Bres = convert(Real, func(low_args...))..convert(Real, func(high_args...))
    return Interval{Float64}(Bres)
end

"""
    @extend_monotonic f(dir [,...])

Define `2^N-1` new funcitons (where `N` is the number of given arguments of `f`)
extending the elementwise monotonic function `f` for intervals computations for
any combination of `Real` and `Interval{Float64}` arguments.

Arguments `dir` of `f` are the direction of the monoticity for that argument,
either `:up` or `:down`.
"""
macro extend_monotonic(expr)
    expr_dict = matchex(:(func(dirs...)), expr, phs=[:func, :dirs])
    expr_dict = get(expr_dict)

    N = length(expr_dict[:dirs])
    func = expr_dict[:func]

    args = [Symbol("X$i") for i in 1:N]
    argtype = Interval{Float64}

    typed_args = [subs(:(A::T), A=arg, T=argtype) for arg in args]
    info_args = [subs(:((A, M)), A=arg, M=monot) for (arg, monot) in zip(args, expr_dict[:dirs])]

    Ifunc_expr = quote
        function $(esc(func))($(typed_args...))
            monotonic_extender($(esc(func)), $(info_args...))
        end
    end


    typed_reals = [subs(:(A::Real), A=arg) for arg in args]
    intervalled_args = [subs(:(Interval(A)), A=arg) for arg in args]
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
