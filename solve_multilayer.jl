push!(LOAD_PATH,pwd())

using PyPlot
using StaticArrays
using ValidatedNumerics

using GeneratingFunctions

const UNITINTERVAL = 0..1

low(interval) = interval.lo
high(interval) = interval.hi

struct Layer
    dist::Symbol
    param::Real
    g0::Function
    g1::Function
    dg0::Function
    dg1::Function

    function Layer(dist::Symbol, param::Real)
        g0 = z -> GFUNC[dist][:g0](z, param)
        g1 = z -> GFUNC[dist][:g1](z, param)
        dg0 = z -> GFUNC[dist][:dg0](z, param)
        dg1 = z -> GFUNC[dist][:dg1](z, param)

        return new(dist, param, g0, g1)
    end
end

"""
    Generate residual function to find the ``u_i`` for a given vector of layers.
"""
function res_gen(uu, layers)
    L = length(layers)
    res = Vector{Any}(L)

    for (j, layer) in enumerate(layers)
        p = 1
        for (i, layer) in enumerate(layers)
            if i == j
                p *= 1 - layer.g1(uu[i])
            else
                p *= 1 - layer.g0(uu[i])
            end
        end

        res[j] = 1 - uu[j] - p
    end

    return SVector{L}(res)
end


function find_u(layers::Vector{Layer}, maxiter=1000)
    L = length(layers)
    X = IntervalBox(UNITINTERVAL, L)
    rts = roots(uu -> res_gen(uu, layers), X)

    affined_rts = []

    for rt in rts
        no_root = false

        if rt.status != :unique
            prev_interval = Interval(-1)
            for iter in 1:maxiter
                prev_interval = rt.interval
                new_rts = roots(uu -> res_gen(uu, layers), rt.interval)

                if length(new_rts) == 0
                    # No solution in the interval
                    no_root = true
                    break
                elseif length(new_rts) == 1
                    rt = new_rts[1]
                    if rt.status == :unique
                        # Solution is unique
                        break
                    elseif prev_interval == rt.interval
                        # Continuing iterating would give nothing
                        break
                    end
                else
                    error("Multiple intervals created. Implementation is missing.")
                end
            end
        end

        if !no_root
            push!(affined_rts, rt)
        end

    end
    return affined_rts
end

function find_S(layers::Vector{Layer}, maxiter=1000)
    L = length(layers)
    uu = find_u(layers)
    SS = []

    for rt in uu
        s = 1
        for (layer, u) in zip(layers, rt.interval)
            s *= 1 - layer.g0(u)
        end
        push!(SS, (s, rt.status))
    end

    return SS
end

function boundary_res_gen(uu_cc, layer_types::Vector{Symbol})
    L = length(layers)
    res = Vector{Any}(length(uu_cc))
    uu = uu_cc[1:L]

    for (j, layer) in enumerate(layers)
        p = 1
        for (i, layer) in enumerate(layers)
            if i == j
                p *= 1 - layer.g1(uu[i])
            else
                p *= 1 - layer.g0(uu[i])
            end
        end

        res[j] = 1 - uu[j] - p
    end

    return SVector{L}(res)
end

function find_boundary(layer_params::Vector{Tuple{Symbol, Interval}}, maxiter=1000)
end

info("Initialization done.")

H = []
M = []
L = []
cc = 1.5:0.01:3
for c in cc
    layers = [Layer(:geometric, c), Layer(:geometric, c)]
    SS = find_S(layers)
    sols = [S[1] for S in SS]
    push!(H, maximum(sols))
    push!(L, minimum(sols))
    if length(sols) > 1
        push!(M, sort(sols)[2])
    else
        push!(M, sols[1])
    end
end

plot(cc, mid.(H), marker=".", linestyle="-", linewidth=0.5)
