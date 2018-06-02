"""
    erdos_renyi(n::Int, c::Real)

Generate an Erdos-Rényi graph with `n` vertices and mean degree `c`.
"""
function erdos_renyi(n::Int, c::Real)
    ne = round(Int, n*c/2) # Number of edges
    v = 1:n
    g = Graph(n)

    for i in 1:ne
        add_edge!(g, Edge(rand(v), rand(v)))
    end

    return g
end

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

    stubs = Vector{Int}(total)

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

abstract type GraphGenerator end

struct ErdosRenyiGraph <: GraphGenerator end
ErdosRenyiGraph(n, c) = erdos_renyi(n, c)

struct ScaleFreeGraph <: GraphGenerator end
ScaleFreeGraph(n, α) = configuration_model(plrand(α, n))

struct GeometricGraph <: GraphGenerator end
GeometricGraph(n, c) = configuration_model(rand(Geometric(1/c), n) + 1)

struct MultiGraph <: GraphGenerator
    layers::Vector{GraphGenerator}
end

MultiGraph(n, parameters) = [layer(n, p) for (layer, p) in zip(mg.layers, parameters)]
