import JSON

mutable struct Simulation{D <: Union{Symbol, Vector{Symbol}}, P, R}
    distribution::D
    n::Int
    parameters::Vector{P}
    repeat::Int
    results::Vector{R}
end

function JSON.lower(sim::Simulation)
    fields = fieldnames(Simulation)
    dict = Dict()
    for field in fields
        dict[String(field)] = getfield(sim, field)
    end
    return dict
end

function save(file, sim::Simulation, replace=false)
    content = []
    if !replace
        if isfile(file)
            print(readlines(file))
            line = readlines(file)[1] # File should always be one line
            for s in JSON.parse(line)
                push!(content, s)
            end
        end
    end
    push!(content, sim)
    write(file, JSON.json(content))
end
