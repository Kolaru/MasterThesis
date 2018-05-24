module Graphs

# TODO Write proper documentation for the file #doc
"""
    Graph{T}

Graph data structure. Have similar API and code as the `LightGraphs.jl`
package but allows for multi-edges.

Was not build on top of `LightGraphs` since the API changed in the last
version which requires Julia 0.7 with which I didn't want to wrestle.
"""
struct Graph{T}
    ne::T
    adjlist::Vector{Vector{T}}
end

struct Edge{T}
    src::T
    dst::T
end

eltype(::Graph{T}) where T = T
edgetype(::Graph{T}) where T = Edge{T}

nv(g::Graph) = length(g.adjlist)
ne(g::Graph) = g.ne
vertices(g::Graph) = one(T):nv(g)

copy(g::Graph) = Graph(g.ne, deepcopy(g.adjlist))

neighbors(g::Graph, v::Integer) = g.adjlist[v]

function add_edge!(g::Graph{T}, edge::Edge{T}) where T
    s = edge.src
    d = edge.dst
    insert_sorted!(g.adjlist[s], d)
    insert_sorted!(g.adjlist[d], s)
end

function rem_edge!(g::Grapht{T}, edge::Edge{T}) where T
    s = edge.src
    d = edge.dst
    remove_sorted!(g.adjlist[s], d)
    remove_sorted!(g.adjlist[d], s)
end

function add_vertex!(g::Graph{T}) where T
    push!(g.adjlist, T[])
end

function rem_vertex!(g::Graph, v::Integer)
    last = nv(g)

    # Remove all edges of v
    for w in neighbors(g, v)
        rem_edge!(g, Edge(v, w))
    end

    # Remove all edges of the last vertex and put them at `v`
    for w in neighbors(g, last)
        rem_edge!(g, Edge(last, w))
        add_edge!(g, Edge(v, w))
    end

    # Delete the last vertex
    pop!(g.adjlist)
end

function in_sorted(vec::Vector{T}, x::T) where T
    searchsortedfirst(vec, x) > length(vec)
end

#= The next function is inspired by
    https://stackoverflow.com/questions/25678112/insert-item-into-a-sorted-list-with-julia-with-and-without-duplicates
=#
function insert_sorted!(vec::Vector{T}, x::T) where T
    insert!(vec, searchsortedfirst(vec, x), x)
end

function replace_sorted!(vec::Vector{T}, xold::T, xnew::T}) where T
    sp = splice!(vec, searchsorted(vec, xold))
    sf = searchsortedfirst(vec, xnew)
    for i in 1:length(sp)
        insert!(v, sf, xnew)
    end
end

function remove_sorted!(vec::Vector{T}, x::T) where T
    splice!(vec, searchsortedfirst(vec, x))
end

function removeall_sorted!(vec::Vector{T}, x::T) where T
    splice!(vec, searchsorted(vec, x))
end

end
