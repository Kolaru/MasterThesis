import StatsBase: fit, Histogram

function poisson(c, k)
    p = 1.
    for j in 1:k
        p *= c/j
    end
    return p*exp(-c)
end

function poisson_no_orphan(c, k)
    k == 0 && return 0.0
    return poisson(c, k)/(1 - exp(-c))
end

function gccpoisson(c, k)
    u = 0.
    for _ in 1:1000
        u = exp(c*(u - 1))
    end
    return poisson(c, k)*(1 - u^k)
end

function mu(r, z)
    N = length(r)
    res = 0.
    Z = 0.
    pik = 1.
    for k in N:-1:1
        pik = r[k]/(1 - z^k)
        Z += k*pik
        res += k*pik*z^(k-1)
    end
    return res/Z
end

function find_global_dist(r, u=0.0, tol=1e-14)
    r[1] == 0 && return 0.0, r
    u = r[1]
    k = 1
    while true
        ustart = u
        u = mu(r, ustart)
        abs(u - ustart)/u < tol && break

        k += 1
        k > 10000 && error("Convergence  not achieved")
    end
    pks = r./(1 .- u.^(1:length(r)))
    return u, pks/sum(pks)
end

function compare_real_with_generated(real_net_name)
    g = RealGraph(real_net_name)
    deg = degrees(g)
    hist = fit(Histogram, deg, nbins=maximum(deg), closed=:left)
    hist = normalize(hist)
    rk = hist.weights

    u, pk = find_global_dist(rk, 0.0, 1e-14)
    ks = 1:length(rk)
    S = 1 - sum(pk .* u.^ks)
    N = 10000000

    gshuffle = configuration_model(deg)
    rkshuffle = gcc_degree_dist(gshuffle)
    rkshuffle = rkshuffle[2:end]  # Exclude degree 0

    Nk = round.(Int, pk*N)
    ggen = EmpiricalGraph(Nk, lowest_degree=1)
    rkgen = gcc_degree_dist(ggen)
    rkgen = rkgen[2:end]  # Exclude degree 0

    println(sum(rk))
    println(sum(pk))
    println(sum(rkgen))
    println(sum(rkshuffle))

    print(gcc_size(gshuffle)/nv(gshuffle))

    nmin = min(length(rk), length(rkgen), length(rkshuffle))

    data = Dict("N" => N,
                "rk" => rk,
                "rkgen" => rkgen,
                "rkshuffle" => rkshuffle
                )
    open("Plot generation/connected_network_generator/$real_net_name.json", "w") do file
        write(file, JSON.json(data))
    end

    scatter(ks[1:nmin], abs.(rkgen[1:nmin] ./ rk[1:nmin]))
    scatter!(ks[1:nmin], abs.(rkshuffle[1:nmin] ./ rk[1:nmin]))

end


function gcc_degree_dist(g)
    components = connected_components(g)
    gcc_ids = components[indmax(length.(components))]
    gcc = subgraph(g, gcc_ids)
    rks = degree_dist(gcc)
    return unshift!(rks, 0.0)
end

function degree_dist(g)
    deg = degrees(g)
    hist = fit(Histogram, deg, nbins=maximum(deg), closed=:left)
    hist = normalize(hist)
    return hist.weights
end

function bins(hist::Histogram)
    edges = hist.edges[1]
    hist.closed == :left && return edges[1:end-1]
    return edges[2:end]
end

function test_gcc_with_poisson()
    K = 1000
    ks = 1:K
    cs = 1.1:0.1:1.4

    all_data = []

    for c in cs
        u = 0.0
        for _ in 1:1000
            u = exp(c*(u - 1))
        end

        rk = gccpoisson.(c, ks)
        rk /= sum(rk)
        urecon, pk = find_global_dist(rk)

        data = Dict("c" => c,
                    "u" => u,
                    "rk" => rk,
                    "urecon" => urecon,
                    "pk" => pk)

        push!(all_data, data)
    end

    open("Plot generation/connected_network_generator/ER_verif.json", "w") do file
        write(file, JSON.json(all_data))
    end
end
