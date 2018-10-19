export Graph, Edge

# TODO Write proper documentation for the file #doc
"""
    Graph{T}

Graph data structure. Have similar API and code as the `LightGraphs.jl`
package but allows for multi-edges.

Was not build on top of `LightGraphs` since the API changed in the last
version which requires Julia 0.7 with which I didn't want to wrestle.
"""
mutable struct Graph{T <: Integer}
    ne::T
    adjlist::Vector{Vector{T}}
end

Graph(n::T) where {T <: Integer} = Graph(0, [Vector{Int}() for _ in 1:n])

struct Edge{T}
    src::T
    dst::T
end

export src, dst, weight

src(e::Edge) = e.src
dst(e::Edge) = e.dst
weight(e::Edge) = 1

export eltype, edgetype, nv, ne, vertices, degrees, copy, neighbors, edges,
       adjacency_matrix, laplacian_matrix

eltype(::Graph{T}) where T = T
edgetype(::Graph{T}) where T = Edge{T}

nv(g::Graph) = length(g.adjlist)
ne(g::Graph) = g.ne
vertices(g::Graph{T}) where T = one(T):nv(g)
degrees(g::Graph) = length.(g.adjlist)
copy(g::Graph) = Graph(nv(g), deepcopy(g.adjlist))
neighbors(g::Graph, v::Integer) = g.adjlist[v]

edges(g::Graph) = Channel() do channel
    for v in vertices(g)
        for n in neighbors(g, v)
            n > v && push!(channel, Edge(v, n))
        end
    end
end

function adjacency_matrix(g::Graph)
    mat = spzeros(Int, nv(g), nv(g))
    for e in edges(g)
        mat[src(e), dst(e)] += 1
        mat[dst(e), src(e)] += 1
    end
    return mat
end

function laplacian_matrix(g::Graph)
    deg = degrees(g)
    lapmat = spzeros(Int, nv(g), nv(g))
    for i in vertices(g)
        lapmat[i, i] = deg[i]
    end
    return lapmat - adjacency_matrix(g)
end

export add_edge!, add_vertex!, grow!, rem_edge!, rem_vertex!

function add_edge!(g::Graph{T}, edge::Edge{T}) where T
    s = edge.src
    d = edge.dst
    insert_sorted!(g.adjlist[s], d)
    s != d && insert_sorted!(g.adjlist[d], s)  # Add only one end for self loop
    g.ne += 1
end

function rem_edge!(g::Graph{T}, edge::Edge{T}) where T
    s = edge.src
    d = edge.dst
    remove_sorted!(g.adjlist[s], d)
    s != d && remove_sorted!(g.adjlist[d], s)  # Remove only one end for self loop
    g.ne -= 1
end

function add_vertex!(g::Graph{T}) where T
    push!(g.adjlist, T[])
end

function grow!(g::Graph{T}, n::Integer) where T
    n <= nv(g) && return
    for k in (nv(g)+1):n
        add_vertex!(g)
    end
end

function rem_vertex!(g::Graph, v::Integer)
    last = nv(g)

    # Remove all edges of v
    for w in deepcopy(neighbors(g, v))
        rem_edge!(g, Edge(v, w))
    end

    processing_self_edge = false
    # Remove all edges of the last vertex and put them at `v`
    for w in deepcopy(neighbors(g, last))
        if w == last
            if !processing_self_edge
                processing_self_edge = true
                rem_edge!(g, Edge(last, w))
                add_edge!(g, Edge(v, v))
            else
                processing_self_edge = false
            end
        else
            rem_edge!(g, Edge(last, w))
            add_edge!(g, Edge(v, w))
        end
    end

    # Delete the last vertex
    pop!(g.adjlist)
end

export subgraph, subgraph!

"""
    subgraph(g::Graph, indices::Vector{Int})

Return the subgraph where only the vertices with indices in `indices` are present.

Vertex `i` in the subgraph has index `indices[i]` in the parent graph. To allow
this behavior, `indices` must be a sorted list.
"""
function subgraph(g::Graph, indices::Vector{Int})
    subg = copy(g)
    subgraph!(subg, indices)
    return subg
end

function subgraph!(g::Graph, indices::Vector{Int})
    n = nv(g)
    for i in setdiff(n:-1:1, indices)
        rem_vertex!(g, i)
    end
end
