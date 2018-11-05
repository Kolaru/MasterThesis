module MonotonicExtension

using Espresso
using IntervalArithmetic

import Base: widen

export AtomicExtension, Extension, Singularity
export @monotone

"""
    widen(X::Interval, relerr)

Widen an interval in order to take in account a relative error of `relerr` on
the computation of its bounds.
"""
function widen(X::Interval, relerr)
    relerr = Interval(relerr)  # Use interval to guarantee computation
    low_mod = X.lo > 0 ? 1 - relerr : 1 + relerr
    high_mod = X.hi > 0 ? 1 + relerr : 1 - relerr

    low_bound = (low_mod*X.lo).lo
    high_bound = (high_mod*X.hi).hi

    return Interval(low_bound, high_bound)
end

# TODO: special case where X totally outsie of clampto
"""
    clamp(X::Interval, clampto::Interval)

Clamp an interval `X` to an interval `clampto`.
"""
function clamp(X::Interval, clampto::Interval)
    clampto == -∞..∞ && return X
    low_bound = max(X.lo, clampto.lo)
    high_bound = min(X.hi, clampto.hi)
    return Interval(low_bound, high_bound)
end

struct Singularity{T, I <: Integer}
    value::T
    index::I
    in_domain::Bool
end

function cut(dom::IntervalBox, sing::Singularity)
    X = dom[sing.index]

    sing.value ∉ X && return dom

    domleft = collect(dom)
    domright = copy(domleft)

    if sing.in_domain
        leftval = rightval = sing.value
    else
        leftval = prevfloat(sing.value)
        rightval = nextfloat(sing.value)
    end

    domleft[sing.index] = Interval(X.lo, leftval)
    domright[sing.index] = Interval(rightval, X.hi)

    return IntervalBox(domleft), IntervalBox(domright)
end

"""
    AtomicExtension{F <: Function, N}

Representation of an extension for an argument wise monotonic function over a
single domain.

# Fields
    - `f`: The function to extend.
    - `domain`: The domain over which `f` is monotonic.
    - `monots`: Store if `f` is increasing or decreasing for each of its arguments.
    - `clampto`: Result of the function are clamp to this interval.
    - `relerr`: Relative error on the computation of `f`.
"""
struct AtomicExtension{F <: Function, N}
    f::F
    domain::IntervalBox{N, Float64}
    monots::Vector{Symbol}
    clampto::Interval{Float64}
    relerr::Float64
end

function AtomicExtension(f, dom, clampto=-Inf..Inf, relerr=0.)
    monots = infer_monotonicity(f, dom)
    AtomicExtension(f, dom, monots, clampto, relerr)
end

function (ext::AtomicExtension{F, N})(Xs)  where {F, N}
    low_args, high_args = get_args(ext, Xs)

    low_bound = prevfloat(ext.f(low_args...))
    high_bound = nextfloat(ext.f(high_args...))

    # Conversion to Real is performed as some special function always return Complex
    try
        low_bound = convert(Real, low_bound)
        high_bound = convert(Real, high_bound)
    catch err
        if isa(err, InexactError)
            warn("Bounds computed by montone extension are not Real.")
            warn("  Low bound: $low_bound")
            warn("  High bound: $high_bound")
        end
        rethrow(err)
    end

    res = Interval(low_bound, high_bound)
    res = widen(res, ext.relerr)

    return clamp(res, ext.clampto)
end

domain_contains(ext::AtomicExtension, x::Vector{T}) where T = all(x .∈ ext.domain)
domain_contains(ext::AtomicExtension, x) = x ∈ ext.domain
domain_contains(ext::AtomicExtension, X::Region) = X ⊆ ext.domain

"""
    Extension{F <: Function, N}

Extension of a function that is monotonic over an union of domain (the monotonicitiy
may change from one domain to the others).
"""
struct Extension{F <: Function, N}
    atomic_extensions::Vector{AtomicExtension{F, N}}
end

function Extension(f, dom, clampto, relerr, ::Nothing)
    Extension([AtomicExtension(f, dom, clampto, relerr)])
end

function Extension(f, dom, clampto, relerr, sing::Singularity)
    Extension(f, dom, clampto, relerr, [sing])
end

function Extension(f, dom::IntervalBox{N, T}, clampto, relerr,
                   singularities::Vector{Singularity{T, I}}) where {N, T, I}
    #
    domains = [dom]
    for sing in singularities
        new_domains = IntervalBox{N, T}[]
        for d in domains
            append!(new_domains, cut(d, sing))
        end
        domains = new_domains
    end

    return Extension([AtomicExtension(f, d, clampto, relerr) for d in domains])
end

function (extension::Extension{F, N})(Xs::Vararg{Interval{Float64}, N})  where {F, N}
    extind = findfirst(ext -> domain_contains(ext, IntervalBox(Xs)), extension.atomic_extensions)

    # TODO: Proper error
    # TODO: Bisect Xs if the singularities are part of the domain
    extind == nothing && error("No subdomain containing $Xs, domains are $([ext.domain for ext in extension.atomic_extensions])")

    return extension.atomic_extensions[extind](Xs)
end
"""
    get_args(ext::AtomicExtension, Xs)

Return two list the arguments that respectively give low and high bounds in each
of the interval in the list `Xs` for the extension `ext`.
"""
function get_args(ext::AtomicExtension{F, N}, Xs) where {F, N}
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

"""
    infer_monotonicity(func, domain)

For each argument of `func` find out if it is increasing or decreasing over the
domain `domain`.
"""
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

        # TODO: Add error if one of the evaluation is NaN

        if func(low_args...) < func(high_args...)
            monots[i] = :increasing
        else
            monots[i] = :decreasing
        end
    end

    return monots
end


"""
    @monotone f Domain(d1, d2, ...) clampto=-Inf..Inf relerr=0. singularities=nothing

Define `2^N-1` new funcitons (where `N` is the number of arguments of `f`)
extending the argumentwise monotonic function `f` for intervals computations for
any combination of `Any` and `Interval{Float64}` arguments.
"""
macro monotone(f, dom, clampto=-Inf..Inf, relerr=0., singularities=nothing)
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
        ext = Extension($func, dom, $clampto, $relerr, $singularities)
        FUNCDEFS...
    end
    return subs(expr, DOMS=doms, FUNCDEFS=funcdefs)
end

end
