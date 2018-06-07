import Distributions: Geometric
import JSON
import PowerLawDistribution: plrand
import StatsBase: sample
using ValidatedNumerics

if !isdefined(:first_run)
    first_run = false
    include("graphs.jl")
    include("network_generation.jl")
    include("simulation.jl")
    include("connectivity.jl")
    include("generating_functions.jl")
end

function run_all()
    n = 1000000
    repeat = 10

    cc = 1.1:0.1:2
    sim = GCCSimulation(GeometricGraph, n, cc, repeat)
    run_simulation!(sim)
    save("Geometric.json", sim, true)
end

# @time run_all()
