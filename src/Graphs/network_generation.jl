export ErdosRenyiGraph, GeometricGraph, ScaleFreeGraph, RealGraph,
    SaturatedScaleFreeGraph, MultiGraph, GraphType, ConnectedGraph

"""
    configuration_model(degrees)

Generate a network given a set of degrees using the configuration model.
"""
function configuration_model(degrees)
    total = sum(degrees)
    # If the total degree is uneven, add a vertex with degree 1
    if total % 2 != 0
        total += 1
        push!(degrees, 1)
    end
    stubs = Vector{Int}(undef, total)

    # Create a list of stubs where vertice v appears deg(v) times
    s = 1
    for (v, d) = enumerate(degrees)
        stubs[s:s+d-1] = fill(v, d)
        s += d
    end
    shuffle!(stubs)

    # Connect the stubs two by two and put the edges in the graph
    g = Graph(length(degrees))
    for i in 1:2:total
        add_edge!(g, Edge(stubs[i], stubs[i+1]))
    end

    return g
end

abstract type GraphType end

struct ConnectedGraph <: GraphType end

struct ErdosRenyiGraph <: GraphType end
function ErdosRenyiGraph(n::Int, c::Real)
    ne = round(Int, n*c/2) # Number of edges
    v = 1:n
    g = Graph(n)

    for i in 1:ne
        add_edge!(g, Edge(rand(v), rand(v)))
    end

    return g
end

struct ScaleFreeGraph <: GraphType end
function ScaleFreeGraph(n, α)
    degs = plrand(α, n)
    maximum(degs) > n && return ConnectedGraph
    return configuration_model(degs)
end

struct GeometricGraph <: GraphType end
GeometricGraph(n, c) = configuration_model(rand(Geometric(1/c), n) .+ 1)

struct SaturatedScaleFreeGraph <: GraphType end

struct EmpiricalGraph <: GraphType end
# EmpiricalGraph(n, pk) = configuration_model(sample(1:length(pk), weights(pk), n))
function EmpiricalGraph(Nk::Vector{Int} ; lowest_degree=0)
    deg = vcat([fill(k - 1 + lowest_degree, Nk[k]) for k in 1:length(Nk)]...)
    return configuration_model(deg)
end

struct RealGraph <: GraphType end
function RealGraph(name)
    path = "Data/real-networks/$name/out.$name"
    g = Graph(1)

    edges = readdlm(path, Int, comments=true, comment_char='%', use_mmap=true)
    for k in 1:size(edges, 1)
        v1, v2 = edges[k, :]
        grow!(g, max(v1, v2))
        add_edge!(g, Edge(v1, v2))
    end
    return g
end

struct MultiGraph <: GraphType end

MultiGraph(n, layers, parameters) = [layer(n, p) for (layer, p) in zip(layers, parameters)]
