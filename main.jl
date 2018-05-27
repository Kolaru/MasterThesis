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
