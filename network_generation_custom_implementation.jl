using DataStructures
using Distributions
using NLsolve
using ProgressMeter
using PyPlot
using SpecialFunctions
using StatsBase



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
