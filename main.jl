import LightGraphs: connected_components,
                    erdos_renyi,
                    neighbors,
                    nv,
                    rem_vertex!

import PowerLawDistribution: plrand

# TODO ensure that the degrees sequences are graphical
# TODO use LightGraphs scale free graph generator #decision
const GRAPHS = Dict(
    :poisson => (n, c) -> erdos_renyi(n, c/n),
    :powerlaw => (n, α) -> begin
        degrees = plrand(α, n)
        if isodd(sum(degrees))
            degrees[1] += 1
        end
        return random_configuration_model(n, degrees, check_graphical=true)
    end,
    :geometric => (n, c) -> error("Not implemented")
)

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

function connected_components(g::Graph, present_indices::Vector{Int})
    n = nv(g)
    subg = subgraph(g, present_indices)
    sub_components = connected_components(subg)
    return [present_indices[comp] for comp in sub_components]
end

function gcc_on_range(degree_dist, n, cc, repeat)
    gcc_sizes = []
    for c in cc
        sizes = Vector{Int}()
        for i in 1:repeat
            net = GRAPHS[degree_dist](n, c)
            comps = connected_components(net)
            push!(sizes, maximum(length.(comps)))
        end
        push!(gcc_sizes, mean(sizes))
    end
    return gcc_sizes/n
end

"""
    extended_neighborhood(g::Graph, start_vertex::Int, present_indices)

Find all vertices connected to vertex `start` in the subgraph of `g` containing
the vertices with indices `present_indices`.
"""
function extended_neighborhood(g::Graph, start_vertex::Int, present_indices)
    n = nv(g)
    is_present = falses(n)
    is_present[present_indices] = true
    queued = falses(n)
    queued[start_vertex] = true

    neigs = Int[start_vertex]
    connected = Int[]

    while !isempty(neigs)
        v = pop!(neigs)
        if is_present[v]
            push!(connected, v)
            for v2 in neighbors(g, v)
                if is_present[v2] && !queued[v2]
                    push!(neigs, v2)
                    queued[v2] = true
                end
            end
        end
    end

    return connected
end


"""
    viable_components_size(multi_net::Vector{Graph})

Find viable components size in the multi layer network `multi_net`.

Algorithm adapted from Baxter 2012.
"""
function viable_components_size(multiplex_network::Vector{Graph{Int}}, fraction_thres=0.001)
    L = length(multiplex_network)  # Number of layers
    n = nv(multiplex_network[1])  # Number of vertices
    n_thres = n*fraction_thres

    components_size = Int[]
    multi_comp = connected_components.(multiplex_network)
    sort!.(multi_comp, by=length)
    multi_gcc = last.(multi_comp)
    viable_guess = fast_intersect(multi_gcc)
    multi_net = [subgraph(net, viable_guess) for net in multiplex_network]

    println("Size of individual GCC: $(length.(multi_gcc))")
    println("Size of GCC intersection: $(length(viable_guess))")

    to_delete = Int[]
    some_deleted = true
    while some_deleted
        some_deleted = false
        for net in multi_net
            components = connected_components(net)
            for comp in components
                if length(comp) < n_thres
                    some_deleted = true
                    for i in comp
                        rem_vertex!.(multi_net, i)
                    end
                end
            end
        end
    end

    println("Size of clean GCC intersection: $(length(viable_guess))")

    ktest = 0
    n = nv(multi_net[1])
    initial_guess = 1:n
    viable_guess = collect(1:n)
    println(n)
    processed = falses(n)
    previous = Int[]

    for i in initial_guess
        if !processed[i]
            first_pass = true
            viable = viable_guess
            ktest += 1
            while first_pass || length(previous) != length(viable)
                first_pass = false
                previous = viable
                for net in multi_net
                    viable = extended_neighborhood(net, i, viable)
                end
            end
            processed[viable] = true
            viable_guess = setdiff(viable_guess, viable)
            push!(components_size, length(viable))
        end
    end

    println("Number of pass : $ktest")
    if !isempty(components_size)
        println("Final GCC size : $(maximum(components_size))")
    else
        println("Final GCC size : 1")
    end
    println()
    return components_size
end

# Going to sets and then sorting is much faster than intersecting vectors
function fast_intersect(ss)
    sets = Set.(ss)
    inters = intersect(sets...)
    return sort(collect(inters))
end
