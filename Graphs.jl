import Base: copy, eltype

include("sorted_utils.jl")

# TODO Write proper documentation for the file #doc
# TODO Add LightGraph dependancy for integratio nwith GraphPlot.jl
"""
    Graph{T}

Graph data structure. Have similar API and code as the `LightGraphs.jl`
package but allows for multi-edges.

Was not build on top of `LightGraphs` since the API changed in the last
version which requires Julia 0.7 with which I didn't want to wrestle.
"""
struct Graph{T <: Integer}
    ne::T
    adjlist::Vector{Vector{T}}
end

Graph(n::T) where {T <: Integer} = Graph(0, [Vector{Int}() for _ in 1:n])

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

neighbors(g::Graph, v::Integer) = deepcopy(g.adjlist[v])

function add_edge!(g::Graph{T}, edge::Edge{T}) where T
    s = edge.src
    d = edge.dst
    insert_sorted!(g.adjlist[s], d)
    insert_sorted!(g.adjlist[d], s)
end

function rem_edge!(g::Graph{T}, edge::Edge{T}) where T
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
    issorted(indices) || error("Subgraph: Indices list should be sorted.")
    n = nv(g)
    for i in setdiff(n:-1:1, indices)
        rem_vertex!(g, i)
    end
end
