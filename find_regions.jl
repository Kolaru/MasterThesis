using IntervalArithmetic
using IntervalRootFinding
using NLsolve
using StaticArrays

import ForwardDiff: jacobian
import JSON
import JSON: json

using GeneratingFunctions
using Graphs
using Utils

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

    z1func(z1, z2, λ1, λ2) = 1 - (1 - g11(z1, λ1))*(1 - g02(z2, λ2))
    z2func(z1, z2, λ1, λ2) = 1 - (1 - g01(z1, λ1))*(1 - g12(z2, λ2))
    Ufunc(U, C) = SVector(z1func(U..., C...), z2func(U..., C...))

    function detJ(z1, z2, λ1, λ2)
        j11 = 1 - dg11(z1, λ1)*(1 - g02(z2, λ2))
        j12 = - (1 - g11(z1, λ1))*dg02(z2, λ2)
        j21 = - dg01(z1, λ1)*(1 - g12(z2, λ2))
        j11 = 1 - (1 - g01(z1, λ1))*dg12(z2, λ2)
        return det([j11 j12 ; j21 j11])
    end

    function B!(res, z1, z2, λ1, λ2)
        res[1:2] = [z1, z2] - Ufunc([z1, z2], [λ1, λ2])
        res[3] = detJ(z1, z2, λ1, λ2)
    end

    numres = zeros(length(x_parameters))
    for (i, c) in enumerate(x_parameters)
        sol = nlsolve((res, UC) -> B!(res, UC..., c), [0, 0, 2.])
        numres[i] = sol.zero[3]
    end
    return numres
end

struct ParamRegion{N, T <: Real}
    region::IntervalBox{N, T}
    roots::Vector{Root{IntervalBox{N, T}}}
end

contains_trivial(rt::Root{IntervalBox{N, T}}) where {N, T} = ones(N) ∈ rt.interval

function nontrivial_sol(rt::Root{IntervalBox{N, T}}) where {N, T}
    has_sol = (rt.status == :exist || rt.status == :unique)
    return has_sol && !contains_trivial(rt)
end

approx_trivial(rt, ε) = all([Z.lo > 1 - ε for Z in rt.interval])
is_trivial(rt, ε) = isunique(rt) && contains_trivial(rt)

function find_region_slow(ψ, Λ0::IntervalBox{N, T}, δ, ε) where {N, T}
    working = [Λ0]
    I = IntervalBox{N, T}
    R = Root{I}
    trivial = I[]
    nontrivial = I[]
    unkown = I[]
    kiter = 0

    while !isempty(working)
        if kiter % 10 == 0
            print("\rState:   $(length(working)) W    $(length(trivial)) T     $(length(nontrivial)) NT    U $(length(unkown))     ")
        end
        kiter += 1

        Λ = pop!(working)

        if diam(Λ) < δ
            push!(unkown, Λ)
        else
            function fλ(z)
                psi = ψ(z, Λ)
                return SVector(psi[1] - z[1], psi[2] - z[2])
            end

            X = IntervalBox(UNITINTERVAL, N)

            rts = roots(fλ, X, Krawczyk, diam(Λ))

            if any([nontrivial_sol(rt) for rt in rts])
                push!(nontrivial, Λ)
            elseif all([approx_trivial(rt, ε) for rt in rts])
                push!(trivial, Λ)
            else
                Λ1, Λ2 = bisect(Λ)
                push!(working, Λ1)
                push!(working, Λ2)
            end
        end

        if length(working) > 10000
            println()
            warn("Limit on number of interval reached")
            error()
            break
        end
    end

    println()
    println("Final State:   $(length(unkown)) unkown    $(length(trivial)) trivial     $(length(nontrivial)) nontrivial          ")

    return unkown, trivial, nontrivial
end

struct BinaryNode{N, T <: Real}
    id::Int
    Λ::IntervalBox{N, T}
    roots::Vector{Root{IntervalBox{N, T}}}
    verified::Dict{Symbol, Bool}
end

function BinaryNode(k, Λ, rts)
    BinaryNode(k, Λ, rts, Dict(:exist => false,
                               :gap => false,
                               :trivial => false))
end

function BinaryNode(k::Int, Λ::IntervalBox{N, T}) where {N, T}
    BinaryNode(k, Λ, [Root(IntervalBox(UNITINTERVAL, N), :unkown)])
end

function BinaryNode(parentnode::BinaryNode, Λ, side::Symbol)
    BinaryNode(subid(parentnode, side), Λ, parentnode.roots, parentnode.verified)
end



id(bnode::BinaryNode) = bnode.id
parambox(bnode::BinaryNode) = bnode.Λ

function level(bnode::BinaryNode)
    lvl = 0
    lead = id(bnode)
    while lead != 0
        lead = lead >> 1
        lvl += 1
    end
    return lvl
end

Base.parent(bnode::BinaryNode) = id(bnode) >> 1
left_subid(bnode::BinaryNode) = id(bnode) << 1
right_subid(bnode::BinaryNode) = (id(bnode) << 1) + 1

function subid(bnode::BinaryNode, side)
    side == :left && return left_subid(bnode)
    side == :right && return right_subid(bnode)
    error("Subid side $side not valid.")
end

nroots(bnode::BinaryNode) = length(bnode.roots)
IntervalRootFinding.roots(bnode::BinaryNode) = bnode.roots

Base.first(bnode::BinaryNode) = first(roots(bnode))
Base.popfirst!(bnode::BinaryNode) = popfirst!(roots(bnode))
Base.push!(bnode::BinaryNode, rt) = push!(roots(bnode), rt)

struct BinaryTree{N, T}
    nodes::Dict{Int, BinaryNode{N, T}}
    working::Vector{Int}
end

BinaryTree(Λ) = BinaryTree(Dict(1 => BinaryNode(1, Λ)), [1])

has_children(tree::BinaryTree, bnode::BinaryNode) = haskey(tree.nodes, left_subid(bnode))
working(tree::BinaryTree) = tree.working

Base.haskey(tree::BinaryTree) = haskey(tree.nodes)
Base.size(tree::BinaryTree) = maximum(keys(tree.nodes))
Base.first(tree::BinaryTree) = tree[first(working(tree))]
Base.getindex(tree::BinaryTree, k::Int) = tree.nodes[k]
Base.setindex!(tree::BinaryTree, val, k::Int) = (tree.nodes[k] = val)

function Base.push!(tree::BinaryTree, bnode)
    k = id(bnode)
    tree[k] = bnode
    push!(working(tree), k)
end

function process(ψ, X, Λ)
    CX = ψ(X, Λ)
    cap = CX ∩ X
    cap == ∅ && return :empty
    CX ⊆ X && return :exist
    return :unkown
end

"""
    Return the intervals created by the Wien diagram of X and Y.

Intersection of X and Y must not be empty.
"""
function cut(X::Interval, Y::Interval, intstat::Symbol, extstat::Symbol)
    X.lo < Y.lo && return cut(Y, X, intstat, extstat)

    a, b = X.lo, Y.lo
    c, d = extrema([X.hi, Y.hi])

    A = Interval(a, b)
    B = Interval(b, c)
    C = Interval(c, d)

    a == b && b == c && return [(B, :intstat)]
    a == b && return [(A, :extstat), (B, :intstat)]
    b == c && return [(B, :intstat), (C, :extstat)]

    return [(A, :extstat), (B, :intstat), (C, :extstat)]
end

function region_binary_tree(ψ, Λ0::IntervalBox{N, T}, δ, ε) where {N, T}
    tree = BinaryTree(Λ0)

    while !isempty(working(tree))
        bnode = first(tree)
        rt = first(bnode)

        I = interval(rt)
        Λ = parambox(bnode)

        if diam(Λ) < diam(I)
            print("\r Current diam (I): $(diam(I)/δ) δ, nroots = $(length(working(tree))), id = $(id(bnode))            ")

            diam(I) < δ && continue
            popfirst!(bnode)

            X, Y = bisect(I)
            CX, Xstat = process(ψ, X, Λ)
            if Xstat == :empty
                #TODO verifiy if Λ trivial
            elseif Xstat == :exist
                if !contains_trivial(X)
                    #TODO remove Λ from the search
                else
                    push!(bnode, Root(X, Xstat))
                end
            elseif Xstat == :unkown
                push!(bnode, Root(X, Xstat))
            end

            Ystat = process(ψ, Y, Λ)
            Ystat != :empty && push!(bnode, Root(Y, Ystat))
        else
            print("\r Current diam (Λ): $(diam(Λ)/δ) δ, nregions = $(length(working(tree))), id = $(id(bnode))            ")

            popfirst!(working(tree))

            diam(Λ) < δ && continue

            Λ1, Λ2 = bisect(Λ)

            push!(tree, BinaryNode(left_subid(bnode), Λ1, roots(bnode)))
            push!(tree, BinaryNode(right_subid(bnode), Λ2, roots(bnode)))
        end
    end

    return tree
end

function find_region_in_two_layers(dist1, dist2, Λ0, δ=5e-3, ε=1e-3)
    g01 = g0(dist1)
    g11 = g1(dist1)

    g02 = g0(dist2)
    g12 = g1(dist2)

    ψ1(z1, z2, λ1, λ2) = 1. - (1. - g11(z1, λ1))*(1. - g02(z2, λ2))
    ψ2(z1, z2, λ1, λ2) = 1. - (1. - g01(z1, λ1))*(1. - g12(z2, λ2))
    ψ(z, λ) = IntervalBox(ψ1(z..., λ...), ψ2(z..., λ...))

    tree = region_binary_tree(ψ, Λ0, δ, ε)

    println(size(tree))
    error()

    data = Dict("unkown" => U,
                "trivial" => T,
                "nontrivial" => NT)

    open("Plot generation/critical_region/$(dist1)_$dist2.json", "w") do file
        JSON.print(file, data)
    end

    return U, T, NT
end

function find_region_in_monolayer(dist, initial_region)
    Ufunc(u, C) = g1(dist, u, C[1], C[2])

    return find_region(Ufunc, initial_region, Val{1})
end
