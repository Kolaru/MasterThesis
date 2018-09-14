import json
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

def plot_real_net_res(name):
    with open(head + "{}.json".format(name)) as file:
        data = json.load(file)

    rk = np.asarray(data["rk"])
    rkshuffle = np.asarray(data["rkshuffle"])
    rkgen = np.asarray(data["rkgen"])

    K = np.min([len(rk), len(rkshuffle), len(rkgen)])
    rk = rk[:K]
    rkshuffle = rkshuffle[:K]
    rkgen = rkgen[:K]
    ks = np.arange(1, K+1)

    fig, ax = plt.subplots(figsize=(4, 3))

    pgen, = ax.plot(ks, rkgen / rk, "o", label="Our algorithm")
    pshuffle, = ax.plot(ks, rkshuffle / rk, "^", label="Edges shuffling")
    ax.set_xlabel("$k$")
    ax.set_ylabel("$r^{gen}_k/r_k$")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(handles=(pgen, pshuffle))
    fig.tight_layout()
    fig.savefig("LaTeX/Degree distribution in GCC/real_{}.pdf".format(name))


plt.rc("text", usetex=True)
plt.rc("font", size=14)

plot_ER = False
plot_road_net = True
plot_power_grid = True

head = "Plot generation/connected_network_generator/"

if plot_ER:
    with open(head + "ER_verif.json") as file:
        all_data = json.load(file)

    fig, ax = plt.subplots()
    K = 20
    ks = np.arange(K)
    smooth_ks = np.linspace(0, K - 1, 1000)

    plots = []
    markers = "o^sdxd"

    for i, data in enumerate(all_data):
        col = "C" + str(i)
        m = markers[i]
        c = data["c"]
        u = data["u"]
        urecon = data["urecon"]
        pk = np.asarray(data["pk"])
        rk = np.asarray(data["rk"])

        S = 1 - (np.exp(c*(u - 1)) - np.exp(-c))/(1 - np.exp(-c))
        p, = ax.plot(ks + 1, rk[ks]/pk[ks], m, color=col, label="$c={}$".format(c), markersize=8)
        ax.plot(smooth_ks + 1, (1 - u**(smooth_ks+1))/S, color=col)

        plots.append(p)

    ax.set_xlabel("$k$")
    ax.set_ylabel("$r_k/p_k$")
    ax.set_xticks([1, 5, 10, 15, 20])
    ax.legend(handles=plots)

if plot_road_net:
    name = "roadNet-CA"
    plot_real_net_res(name)

if plot_power_grid:
    name = "US-power-grid"
    plot_real_net_res(name)

plt.show()
