module MonotonicExtension

using Espresso
using IntervalArithmetic

import Base: widen

export Extension, @monotone

 function widen(X::Interval, relerr)
     relerr = Interval(relerr)  # Use interval to guarantee computation
     low_mod = X.lo > 0 ? 1 - relerr : 1 + relerr
     high_mod = X.hi > 0 ? 1 + relerr : 1 - relerr

     low_bound = (low_mod*X.lo).lo
     high_bound = (high_mod*X.hi).hi

     return Interval(low_bound, high_bound)
 end

function clamp(X::Interval, clampto::Interval)
    clampto == -∞..∞ && return X
    low_bound = max(res.lo, clamp.lo)
    high_bound = min(res.hi, clamp.hi)
    return Interval(low_bound, high_bound)
end

struct Extension{F <: Function, N}
    f::F
    domain::IntervalBox{N, Float64}
    monots::Vector{Symbol}
    clampto::Interval{Float64}
    relerr::Float64
end

function Extension(f, dom, clampto=-Inf..Inf, relerr=0.)
    monots = infer_monotonicity(f, dom)
    Extension(f, dom, monots, clampto, relerr)
end

"""
    get_args(ext::Extension, Xs)

Return two list the arguments that respectively give low and high bounds in each
of the interval in the list `Xs` for the extension `ext`.
"""
function get_args(ext::Extension{F, N}, Xs) where {F, N}
    lows = zeros(N)
    his = zeros(N)
    for k in 1:N
        X = Xs[k]
        if ext.monots[k] == :increasing
            lows[k] = X.lo
            his[k] = X.hi
        else
            lows[k] = X.hi
            his[k] = X.lo
        end
    end
    return lows, his
end

function (ext::Extension{F, N})(Xs::Vararg{Interval{Float64}, N})  where {F, N}
    low_args, high_args = get_args(ext, Xs)

    low_bound = prevfloat(ext.f(low_args...))
    high_bound = nextfloat(ext.f(high_args...))

    # Conversion to Real is performed as some special function always return Complex
    try
        low_bound = convert(Real, low_bound)
        high_bound = convert(Real, high_bound)
    catch err
        if isa(err, InexactError)
            warn("Bounds computed by montone extension are complex.")
            warn("  Low bound: $low_bound")
            warn("  High bound: $high_bound")
        end
        rethrow(err)
    end

    res = Interval(low_bound, high_bound)
    res = widen(res, ext.relerr)

    return clamp(res, ext.clampto)
end

function infer_monotonicity(func, domain)
    N = length(domain)

    mids = mid(domain)
    los = [dom.lo for dom in domain]
    his = [dom.hi for dom in domain]
    monots = Vector{Symbol}(undef, N)

    for i in 1:N
        low_args = collect(mids)
        high_args = collect(mids)
        low_args[i] = los[i]
        high_args[i] = his[i]

        if func(low_args...) < func(high_args...)
            monots[i] = :increasing
        else
            monots[i] = :decreasing
        end
    end

    return monots
end


"""
    @monotone f Domain(d1, d2, ...) clampto=C domain=D

Define `2^N-1` new funcitons (where `N` is the number of arguments of `f`)
extending the argumentwise monotonic function `f` for intervals computations for
any combination of `Any` and `Interval{Float64}` arguments.
"""
macro monotone(f, dom, clampto=-Inf..Inf, relerr=0.)
    func = esc(f)
    matcheddom = matchex(:(Domain(DOMS...)), dom ; phs=[:DOMS])
    if matcheddom != nothing
        doms = [esc(dom) for dom in matcheddom[:DOMS]]
        N = length(doms)
    else
        throw(ArgumentError("Domain for extension must be given as `Domain(d1, d2, ...)`."))
    end

    args = [Symbol("X$i") for i in 1:N]
    interval_args = [subs(:(x::Interval{Float64}), x=arg) for arg in args]
    intervalled_inputs = [subs(:(Interval(x)), x=arg) for arg in args]

    combinations = Iterators.product([[:bare, :interval] for _ in 1:N]...)
    combinations = reshape(collect(combinations), 2^N)

    funcdef = quote
        $func(ARGS...) = ext(INPUTS...)
    end
    funcdefs = []

    for combi in combinations[2:end]
        mixed_args = Vector(undef, N)
        inputs = Vector(undef, N)
        for (k, c) in enumerate(combi)
            if c == :bare
                mixed_args[k] = args[k]
                inputs[k] = intervalled_inputs[k]
            else
                mixed_args[k] = interval_args[k]
                inputs[k] = args[k]
            end
        end
        push!(funcdefs, subs(funcdef, ARGS=mixed_args, INPUTS=inputs))
    end

    expr = quote
        dom = IntervalBox(DOMS...)
        ext = Extension($func, dom, $clampto, $relerr)
        FUNCDEFS...
    end
    return subs(expr, DOMS=doms, FUNCDEFS=funcdefs)
end

end
