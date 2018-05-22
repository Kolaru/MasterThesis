using Distributions
using JSON

import Base.size, Base.getindex, Base.setindex!

mutable struct Node
    id::Int
    deg_id::Int
    neig::Vector{Int}
end

Node(id::Int, deg_id::Int) = Node(id, deg_id, zeros(Int, degrees[deg_id]))
isfull(node::Node) = all(x -> x > 0, node.neig)
degree(node::Node) = length(node.neig)

struct DegreeList <: AbstractArray{Int, 1}
    data::Vector{Int}
    used::Vector{Bool}
end

function DegreeList(data::Vector{Int})
    if isodd(sum(data))
        push!(data, 1)
    end
    DegreeList(sort(data), fill(true, length(data)))
end

size(DL::DegreeList) = size(DL.data)
getindex(DL::DegreeList, i::Int) = DL.data[i]
setindex!(DL::DegreeList, val::Int, i::Int) = (DL.data[i] = val)

all_used(DL::DegreeList) = all(DL.used)
is_available(DL::DegreeList, i::Int) = !DL.used[i]

function promote!(DL, node)
    old_deg_id = node.deg_id
    new_deg_id = findmin(i -> is_available(DL, i) && degree(node) < DL[id], 1:length(DL))
    DL.used[old_deg_id] = false
    DL.used[new_deg_id] = true
    append!(node.neig, zeros(Int, DL[new_deg_id] - DL[old_deg_id]))
end


c = 1.5
dist = Poisson(c)
degrees = DegreeList(rand(dist, 10) + 1)
degrees.available[1] = false
network = [Node(1, 1)]

ktest = 0
while (!all_used(degrees) || !all(isfull, network)) && ktest < 100
    ktest += 1

    if !all_used(degrees)
        if all(isfull, network)
            promote_random_node!(network, degrees)
        else
            add_node!(network, degrees)
        end
    else
        connect_random_nodes!(network, degrees)
    end
end


function promote_random_node!(network, degrees)
    highest_degree = degrees[findlast(degrees.available, true)]
    promotable = [node for node in network if degree(node) < highest_degree]
    promote!(degrees, rand(promotable))
end

function add_node!(network, degrees)

end

function connect_random_nodes!(network, degrees)
end

open("network.json", "w") do f
    write(f, json(network))
end
