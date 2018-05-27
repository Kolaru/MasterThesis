import LightGraphs: connected_components
import PowerLawDistribution: plrand
import StatsBase: sample

include("Graphs.jl")  # TODO Remove leading capital letter
include("network_generation.jl")
include("connectivity.jl")

function gcc_on_range(distribution, n, cc, repeat)
    gcc_sizes = []
    for c in cc
        sizes = Vector{Int}()
        for i in 1:repeat
            net = GRAPHS[distribution](n, c)
            comps = connected_components(net)
            push!(sizes, maximum(length.(comps)))
        end
        push!(gcc_sizes, mean(sizes))
    end
    return gcc_sizes/n
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
