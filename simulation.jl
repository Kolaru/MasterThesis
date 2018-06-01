abstract type Simulation end

mutable struct GCCSimulation{G <: Type{GG} where GG, P, R} <: Simulation
    generator::G
    n::Int
    parameters::Vector{P}
    repeat::Int
    results::Vector{R}
end

GCCSimulation(gen, n, parameters, repeat=1) =
    GCCSimulation(gen, n, collect(parameters), repeat, Vector{typeof(1.0)}())

# Tell JSON to turn the object into a dictionnary.
function JSON.lower(sim::S) where S <: Simulation
    fields = fieldnames(S)
    dict = Dict()
    for field in fields
        dict[String(field)] = getfield(sim, field)
    end
    return dict
end

function save(file, sim::Simulation, replace=false)
    path = "Data/Simulations/$file"
    content = []
    if !replace
        if isfile(path)
            print(readlines(path))
            line = readlines(path)[1] # File should always be one line
            for s in JSON.parse(line)
                push!(content, s)
            end
        end
    end
    push!(content, sim)
    write(path, JSON.json(content))
end

function run_simulation!(sim::GCCSimulation)
    if sim.generator <: MultiGraph
        run_multi_layer_simulation!(sim)
    else
        run_single_layer_simulation!(sim)
    end
end

function run_single_layer_simulation!(sim::GCCSimulation)
    gcc_sizes = typeof(sim.results)()
    for p in sim.parameters
        sizes = Vector{Int}()
        for i in 1:sim.repeat
            net = sim.generator(sim.n, p)
            comps = connected_components(net)
            push!(sizes, maximum(length.(comps)))
        end
        push!(gcc_sizes, mean(sizes))
    end
    append!(sim.results, gcc_sizes/sim.n)
end

function run_multi_layer_simulation!(sim::GCCSimulation)
    viable_sizes = typeof(sim.results)()
    for plist in sim.parameters
        sizes = Vector{Int}()
        for i in 1:sim.repeat
            net = sim.generator(sim.n, plist)
            viables = viable_components_size(net)
            push!(sizes, maximum(length.(viables)))
        end
        push!(viable_sizes, mean(sizes))
    end
    append!(sim.results, viable_sizes/n)
end
