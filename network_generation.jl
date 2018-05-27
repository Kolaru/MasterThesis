"""
    erdos_renyi(n::Int, c::Real)

Generate an Erdos-Rényi graph with `n` vertices and mean degree `c`.
"""
function erdos_renyi(n::Int, c::Real)
    L = round(Int, n*c)
    v = 1:n
    g = Graph(L)

    for i in 1:L
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
    g = Graph(total)

    for i in 1:2:total
        add_edge!(g, Edge(stubs[i], stubs[i]+1))
    end

    return g
end

const GRAPHS = Dict(
    :poisson => erdos_renyi,
    :powerlaw => (n, c) -> configuration_model(plrand(α, n)),
    :geometric => (n, c) -> configuration_model(rand(Geometric(1/c), n - 1) + 1)
)
