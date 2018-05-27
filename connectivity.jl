"""
    fast_instersect(ss)

Efficiently compute the intersection of the vectors contained in `ss`.

The method is faster thant the default, as going to sets and then sorting is
much faster than intersecting vectors directly.
"""
function fast_intersect(ss)
    sets = Set.(ss)
    inters = intersect(sets...)
    return sort(collect(inters))
end

"""
    extended_neighborhood(g::Graph, start_vertex::Int, present_indices)

Find all vertices connected to vertex `start` in the graph `g`.
"""
function extended_neighborhood(g::Graph, start_vertex::Int)
    n = nv(g)
    queued = falses(n)
    queued[start_vertex] = true

    neigs = Int[start_vertex]
    connected = Int[]

    while !isempty(neigs)
        v = pop!(neigs)
        push!(connected, v)
        for v2 in neighbors(g, v)
            if !queued[v2]
                push!(neigs, v2)
                queued[v2] = true
            end
        end
    end

    return connected
end

"""
    connected_components(g::Graph [, present_indices])

Find the connected components in the graph `g` where only the vertices
with indices in `present_indices` are present. By default all vertices are
present.
"""
function connected_components(base_g::Graph)
    n = nv(g)
    processed = falses(n)

    comp = Int[]
    components = Vector{Int}[]

    for i in 1:n
        if !processed[i]
            connected = extended_neighborhood(g, i)
            push!(components, connected)
            for k in connected
                processed[k] = true
            end
        end
    end
    return components
end

function connected_components(g::Graph, present_indices::Vector{Int})
    n = nv(g)
    subg = subgraph(g, present_indices)
    sub_components = connected_components(subg)
    return [present_indices[comp] for comp in sub_components]
end
