import JSON
import PowerLawDistribution: plrand
import StatsBase: sample

include("graphs.jl")
include("network_generation.jl")
include("simulation.jl")
include("connectivity.jl")

cc = 0:0.1:2
sim = GCCSimulation(ErdosRenyiGraph, 10000, cc, 10)
run_simulation!(sim)
save("Erdos_Renyi.json", sim)
