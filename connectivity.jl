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
    extended_neighborhood!(queued::BitArray, g::Graph, start_vertex::Int)

Find all vertices connected to vertex `start` in the graph `g`. `queued` is
a boolean list of telling if a given vertex has been queued.

The syntax allow to reuse the same array, avoiding to reallocate it.
"""
function extended_neighborhood!(queued::BitArray, g::Graph, start_vertex::Int)
    n = nv(g)
    queued[start_vertex] = true
    neigs = Int[start_vertex]
    connected = Int[]

    while !isempty(neigs)
        v = pop!(neigs)
        push!(connected, v)
        for v2 in neighbors(g, v)
            queued[v2] && continue

            push!(neigs, v2)
            queued[v2] = true
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
function connected_components(g::Graph)
    n = nv(g)
    processed = falses(n)
    queued = falses(n)
    components = Vector{Int}[]

    for i in 1:n
        processed[i] && continue

        connected = extended_neighborhood!(queued, g, i)
        push!(components, connected)
        for k in connected
            processed[k] = true
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


function components_size(g::Graph)
    return length.(connected_components(g))
end

function gcc_size(g::Graph)
    return maximum(components_size(g))
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
                length(comp) >= n_thres && continue

                some_deleted = true
                for i in comp
                    rem_vertex!.(multi_net, i)
                end
            end
        end
    end

    println("Size of clean GCC intersection: $(length(viable_guess))")

    ktest = 0
    n = nv(multi_net[1])
    initial_guess = 1:n
    viable_guess = collect(1:n)
    processed = falses(n)
    previous = Int[]

    for i in initial_guess
        processed[i] && continue

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

    println("Number of pass : $ktest")
    if !isempty(components_size)
        println("Final GCC size : $(maximum(components_size))")
    else
        println("Final GCC size : 1")
    end
    println()
    return components_size
end
