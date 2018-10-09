import Distributions: Geometric
import JSON
import StatsBase: sample
using IntervalArithmetic
using IntervalRootFinding
using Plots

import PowerLawDistribution: plrand

if !isdefined(:first_run)
    first_run = false
    include("graphs.jl")
    include("network_generation.jl")
    include("simulation.jl")
    include("connectivity.jl")
    include("src/generating_functions/generating_functions.jl")
    include("find_regions.jl")
    include("connected_network_generator.jl")
end

function simulate_geometric()
    for (n, rep) in [(100, 10000), (1000, 1000), (1000000, 10)]
        info("Geometric simulation with n = $n")
        cc = 1.1:0.1:2
        sim = GCCSimulation(GeometricGraph, n, cc, rep)
        run_simulation!(sim)
        save("Geometric.json", sim)
    end
end

function simulate_ER()
    for (n, rep) in [(100, 10000), (1000, 1000), (1000000, 10)]
        info("ER simulation with n = $n")
        cc = 0.5:0.1:1.5
        sim = GCCSimulation(ErdosRenyiGraph, n, cc, rep)
        run_simulation!(sim)
        save("ER.json", sim)
    end
end

function simulate_scalefree()
    for (n, rep) in [(100, 10000), (1000, 1000), (1000000, 10)]
        info("Scale free simulation with n = $n")
        aa = 2.5:0.2:4.5
        sim = GCCSimulation(ScaleFreeGraph, n, aa, rep)
        run_simulation!(sim)
        save("Scalefree.json", sim)
    end
end
