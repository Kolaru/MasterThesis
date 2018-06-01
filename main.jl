import JSON
import PowerLawDistribution: plrand
import StatsBase: sample

include("graphs.jl")
include("network_generation.jl")
include("simulation.jl")
include("connectivity.jl")

function run_all()
    n = 100000
    repeat = 10

    cc = 0:0.1:2
    sim = GCCSimulation(ErdosRenyiGraph, n, cc, repeat)
    run_simulation!(sim)
    save("Erdos_Renyi.json", sim, true)
end
