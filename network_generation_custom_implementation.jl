using DataStructures
using Distributions
using NLsolve
using ProgressMeter
using PyPlot
using SpecialFunctions
using StatsBase

"""
    Generate an Erdos-Rényi network with `n` vertices and mean degree `c`.
"""
function poisson_network(n::Int, c::Real)
    L = round(Int, n*c)
    v = 1:n
    I = sample(v, L)
    J = sample(v, L)

    return edges_to_adjacency([I J][J .>= I, :], n)
end

"""
    Generate a scale free network with `n` vertices and mean degree `c`.
"""
function powerlaw_network(n::Int, c::Real)
    f! = (F, a) -> plmean_res!(F, a, c)
    res = nlsolve(f!, [2.1])
    α = res.zero[1]

    return edges_to_adjacency(CM(plrand(α, n)), n)
end

function geometric_network(n::Int, c::Real)
    rng = Geometric(1/c)

    return edges_to_adjacency(CM(rand(Geometric(1/c), n - 1) + 1), n)
end

"""
    Create an uncorrelated multi layer network as a list of networks.

    The multi layer network is reis_presented as list of adjacency.
"""
function multi_network(network_types, n, mean_degrees)
    layers = []

    for (net, c) = zip(network_types, mean_degrees)
        push!(layers, net(n, c))
    end

    return layers
end

"""
    Generate network model given a set of degrees using the configuration model.

    Return a list of edges.
"""
function CM(degrees)
    total = sum(degrees)
    # If the total degree is uneven, add a vertex with degree 1
    if total % 2 != 0
        total += 1
        push!(degrees, 1)
    end

    stubs = ones(Int, total)

    # Create a list of stubs where vertice v appears deg(v) times
    s = 1
    for (v, d) = enumerate(degrees)
        stubs[s:s+d-1] = fill(v, d)
        s += d
    end

    shuffle!(stubs)
    half = div(total, 2)
    return [stubs[1:half] stubs[half+1:end]]
end

"""
    Transform a vector of edges to a list of adjacency for a network with n vertices.
"""
function edges_to_adjacency(vertices::Matrix{Int}, n::Int)
    I = vertices[:, 1]
    J = vertices[:, 2]

    if length(I) == 0
        return Vector{Int}[]
    end

    adj = Vector{Vector{Int}}(n)

    for i = 1:n
        adj[i] = Int[]
    end

    for i = 1:length(I)
        push!(adj[I[i]], J[i])
        push!(adj[J[i]], I[i])
    end

    return adj
end

"""
    Find vertices connected to vertex `start` in the network.

    Algortihm from Newmann 2010 (Section 10.3.4).
"""
function connected_vertices(start::Int, network::Vector{Vector{Int}})
    queued = fill(false, length(network))
    queued[start] = true

    neigs = [start]

    connected = Int[]

    while !isempty(neigs)
        v = shift!(neigs)
        push!(connected, v)

        for v2 = network[v]
            if !queued[v2]
                push!(neigs, v2)
                queued[v2] = true
            end
        end
    end

    return connected
end

function connected_vertices(start::Int, network::Vector{Vector{Int}}, indexes, n)
    is_present = fill!(BitArray(n), false)
    @inbounds is_present[indexes] = true
    queued = fill!(BitArray(n), false)
    queued[start] = true

    neigs = Int[start]
    connected = Int[]

    while !isempty(neigs)
        v = shift!(neigs)
        if is_present[v]
            push!(connected, v)
            for v2 in network[v]
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
    Find all components in a network defined as an adjacency list.

    Algortihm from Newmann 2010 (Section 10.3.4).
"""
function find_components(network::Vector{Vector{Int}})
    n = length(network)  # Number of vertices

    processed = fill(false, n)

    comp = Int[]
    components = Vector{Int}[]
    neigs = Queue(Int)

    for i = 1:n
        if !processed[i]
            connected = connected_vertices(i, network)
            push!(components, connected)
            for k = connected
                processed[k] = true
            end
        end
    end
    return components
end

function find_components(network::Vector{Vector{Int}}, indexes, n)
    processed = fill(true, n)
    processed[indexes] = false

    comp = Int[]
    components = Vector{Int}[]
    neigs = Queue(Int)

    for i in indexes
        if !processed[i]
            connected = connected_vertices(i, network, indexes, n)
            push!(components, connected)
            for k = connected
                processed[k] = true
            end
        end
    end
    return components
end

function find_gcc(network::Vector{Vector{Int}})
    components = find_components(network)
    _, maxind = findmax(length.(components))
    return components[maxind]
end


"""
    Find viable components size in a multi layer network.

    Algorithm adapted from Baxter 2012.
"""
function viable_components_size(multi_net)
    n = length(multi_net[1])  # Number of vertices
    L = length(multi_net)  # Number of layers

    components_size = Int[]
    multi_gcc = find_gcc.(multi_net)
    viable_guess = intersect(multi_gcc...)
    processed = fill!(BitArray(n), false)

    last = Int[]
    ktest = 0
    println("GCC intersection size = $(length(viable_guess))")

    to_delete = Int[]
    first_pass = true
    while first_pass || !isempty(to_delete)
        first_pass = false

        filter!(v -> v ∉ to_delete, viable_guess)
        empty!(to_delete)

        for net in multi_net
            components = find_components(net, viable_guess, n)
            for comp in components
                if length(comp) < n*MINIMAL_FRACTION
                    push!(to_delete, comp...)
                end
            end
        end
    end

    initial_guess = copy(viable_guess)
    println("Clean GCC intersection size = $(length(viable_guess))")

    for i in initial_guess
        if !processed[i]
            first_pass = true
            viable = viable_guess
            ktest += 1
            while first_pass || length(last) != length(viable)
                first_pass = false
                last = viable
                for net in multi_net
                    viable = connected_vertices(i, net, viable, n)
                end
            end
            for j in viable
                processed[j] = true
            end
            filter!(v -> v ∉ viable, viable_guess)
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

function viable_components_size_intersections(multi_net)
    n = length(multi_net[1])  # Number of vertices
    L = length(multi_net)

    processed = fill(false, n)
    components_size = Int[]
    multi_components = []

    for network = multi_net
        push!(multi_components, find_components(network))
    end

    lens = [length(comps) for comps = multi_components]
    viable_components = []
    indexes = Base.product(map(n -> 1:n, lens)...)
    short_multi_enum = enumerate(multi_components[2:end])
    for ind = indexes
        viable = multi_components[1][ind[1]]
        for (k, comps) = short_multi_enum
            k += 1
            viable = intersect(viable, comps[ind[k]])
        end

        push!(viable_components, viable)
    end

    return [length(viable) for viable = viable_components]
end


"""
    Find the size of all components in a network defined as an adjacency list.

    Algorithm from Newmann 2010 (Section 10.3.4).
    The algorithm is the same as for `find_components`, but should be faster as
    the list of vertices in each components is not stored.
"""
function find_components_size(network::Vector{Vector{Int}})
    n = length(network)
    processed = fill(false, n)
    queued = fill(false, n)
    c_sizes = Int[]
    c_size = 0
    neigs = Queue(Int)

    for i = 1:n
        if !processed[i]
            c_size = 0
            enqueue!(neigs, i)
            queued[i] = true

            while !isempty(neigs)
                c_size += 1
                v = dequeue!(neigs)
                processed[v] = true

                for v2 = network[v]
                    if !processed[v2] && !queued[v2]
                        enqueue!(neigs, v2)
                        queued[v2] = true
                    end
                end
            end

            push!(c_sizes, c_size)
        end
    end

    return c_sizes
end

function simulate_components_size(network_generator, n, c)
    net = network_generator(n, c)
    return find_components_size(net)
end

function simulate_multi_components_size(generators, n , cs)
    net = multi_network(generators, n, cs)
    vcs = viable_components_size(net)
    return isempty(vcs) ? 1 : maximum(vcs)
end


function simulate_gcc_on_range(network_generator, n, cc, repeat)
    gcc_sizes = []
    for c in cc
        comp_sizes = [simulate_components_size(network_generator, n, c) for i in 1:repeat]
        push!(gcc_sizes, mean(maximum.(comp_sizes)))
    end
    return gcc_sizes/n
end

function simulate_multi_gcc_on_ranges(generators, n, ranges, repeat)
    gcc_sizes = []
    for cs in zip(ranges...)
        comp_sizes = [simulate_multi_components_size(generators, n, cs) for i in 1:repeat]
        push!(gcc_sizes, mean(maximum.(comp_sizes)))
    end

    return gcc_sizes/n
end
