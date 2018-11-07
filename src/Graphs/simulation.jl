export GCCSimulation
export save, run_simulation
export simulate_ER, simulate_geometric, simulate_scalefree

abstract type Simulation end

mutable struct GCCSimulation{G <: Type{GG} where GG, P, R} <: Simulation
    generator::G
    n::Int
    parameters::Vector{P}
    repeat::Int
    sizes::Vector{R}
    stds::Vector{R}
end

mutable struct GVCSimulation{P, R} <: Simulation
    layers::Vector{DataType}
    n::Int
    L::Int
    parameters::Vector{P}
    repeat::Int
    sizes::Vector{R}
    stds::Vector{R}
end

GCCSimulation(gen, n, parameters, repeat=1) =
    GCCSimulation(gen, n, collect(parameters), repeat, Vector{typeof(1.0)}(), Vector{typeof(1.0)}())


GVCSimulation(layers, n, L, parameters, repeat=1) =
    GVCSimulation(layers, n, L, collect(parameters), repeat, Vector{typeof(1.0)}(), Vector{typeof(1.0)}())

# Tell JSON how to turn the object into a dictionnary.
function JSON.lower(sim::S) where S <: Simulation
    fields = fieldnames(S)
    dict = Dict()
    for field in fields
        dict[String(field)] = getfield(sim, field)
    end
    return dict
end

function save(file, sim::Simulation, head="Plot generation/gcc_plots/", replace=false)
    path = head * file
    content = []
    if !replace
        if isfile(path)
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
    gcc_sizes = typeof(sim.sizes)()
    gcc_stds = typeof(sim.sizes)()
    @showprogress 1 for p in sim.parameters
        sizes = Vector{Int}()
        for i in 1:sim.repeat
            net = sim.generator(sim.n, p)
            comps = connected_components(net)
            push!(sizes, maximum(length.(comps)))
        end
        push!(gcc_sizes, mean(sizes))
        push!(gcc_stds, std(sizes))
    end
    append!(sim.sizes, gcc_sizes/sim.n)
    append!(sim.stds, gcc_stds/sim.n)
end

function run_simulation!(sim::GVCSimulation)
    viable_sizes = typeof(sim.sizes)()
    viable_stds = typeof(sim.stds)()
    @showprogress 1 for p in sim.parameters
        sizes = Vector{Int}()
        for i in 1:sim.repeat
            network = MultiGraph(sim.n, sim.layers, fill(p, sim.L))
            filter!(net -> net != ConnectedGraph, network)
            if isempty(network)
                push!(sizes, sim.n)
            else
                viables = viable_components_size(network)
                push!(sizes, maximum(viables))
            end
        end
        push!(viable_sizes, mean(sizes))
        push!(viable_stds, std(sizes))
    end
    append!(sim.sizes, viable_sizes/sim.n)
    append!(sim.stds, viable_stds/sim.n)

    println(sim)
end

function simulate_geometric(L=1)
    if L == 1
        for (n, rep) in [(100, 100000), (1000, 10000), (1000000, 10)]
            @info "Geometric simulation with n = $n"
            cc = 1.1:0.1:2
            sim = GCCSimulation(GeometricGraph, n, cc, rep)
            run_simulation!(sim)
            save("Geometric.json", sim)
        end
    else
        layers = [GeometricGraph for _ in 1:L]
        n = 100000
        rep = 20
        cc = 2:0.3:5
        sim = GVCSimulation(layers, n, L, cc, rep)
        run_simulation!(sim)
        save("GeometricGraph$(L)_sim.json", sim, "Plot generation/single_param_multiplex/", true)
    end
end

function simulate_ER(L=1)
    if L == 1
        for (n, rep) in [(100, 100000), (1000, 10000), (1000000, 10)]
            cc = 0.5:0.1:1.5
            sim = GCCSimulation(ErdosRenyiGraph, n, cc, rep)
            run_simulation!(sim)
            save("ER.json", sim)
        end
    else
        layers = [ErdosRenyiGraph for _ in 1:L]
        n = 100000
        rep = 20
        cc = 2:0.2:4
        sim = GVCSimulation(layers, n, L, cc, rep)
        run_simulation!(sim)
        save("ErdosRenyiGraph$(L)_sim.json", sim, "Plot generation/single_param_multiplex/", true)
    end
end

function simulate_scalefree(L=1)
    if L == 1
        for (n, rep) in [(100, 100000), (1000, 10000), (1000000, 10)]
            aa = 2.5:0.2:4.5
            sim = GCCSimulation(ScaleFreeGraph, n, aa, rep)
            run_simulation!(sim)
            save("Scalefree.json", sim)
        end
    else
        layers = [ScaleFreeGraph for _ in 1:L]
        n = 100000
        rep = 20
        cc = 1.5:0.1:2.5
        sim = GVCSimulation(layers, n, L, cc, rep)
        run_simulation!(sim)
        save("ScaleFreeGraph$(L)_sim.json", sim, "Plot generation/single_param_multiplex/", true)
    end
end
