import json
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

plt.rc("text", usetex=True)
plt.rc("font", size=10, family="palatino linotype")
figsize = (6, 3)

head = "Plot generation/single_param_multiplex/"
markers = "oooooo"

def plot_gvc(ax, name, paramname="c"):
        with open(head + name + ".json") as file:
            data = json.load(file)

        c0 = np.min(data[0]["parameters"])
        c1 = np.max(data[0]["parameters"])

        for k, (d, m) in enumerate(zip(data, markers)):
            ax.errorbar(d["parameters"], d["sizes"], yerr=d["stds"],
                             marker=".",
                             color="C{}".format(k),
                             linewidth=0,
                             capsize=4,
                             elinewidth=2,
                             markeredgewidth=2)
            ax.set_title(label="$n = {}$, $M = {}$".format(d["n"], d["repeat"]))

        ax.set_ylabel("$S$")

def plot_patches(ax, name):
    with open(head + name) as file:
        all_data = json.load(file)

    unkown = []
    unique = []
    exist = []
    xmin, xmax = 1000, 0

    for data in all_data:
        yb, xb = data["bounds"]
        if xb[0] < xmin:
            xmin = xb[0]
        if xb[1] > xmax:
            xmax = xb[1]

        w = xb[1] - xb[0]
        h = yb[1] - yb[0]
        rec = Rectangle([xb[0], yb[0]], w, h)
        if data["status"] == "exist":
            # print("WARNING: data with status 'exist'")
            exist.append(rec)
        elif data["status"] == "unkown":
            unkown.append(rec)
        elif data["status"] == "unique":
            unique.append(rec)

    ax.add_collection(PatchCollection(unique, facecolor="C0"))
    ax.add_collection(PatchCollection(unkown, facecolor="C1"))
    ax.add_collection(PatchCollection(exist, facecolor="C2"))
    return xmin, xmax

def plot_regions(ax, name, n_layers, labelpos, legpos, xlabel):
    xmin, xmax = 1000, 0
    for k, L in enumerate(n_layers):
        fullname = "{}{}.json".format(name, L)
        xa, xb = plot_patches(ax, fullname)

        if xa < xmin:
            xmin = xa
        if xb > xmax:
            xmax = xb

        if labelpos is not None:
            x, y = labelpos[k]
            ax.text(x, y, "$L = {}$".format(L))
                    # bbox=dict(facecolor="white", alpha=0.5, linewidth=0))

    ax.legend((Rectangle((0, 0), 1, 1, facecolor="C0"),
                Rectangle((0, 0), 1, 1, facecolor="C1")),
                ("Unique solution", "Unkown status"), loc=legpos)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("$u$")
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlim(xmin, xmax)

def plot_single_param(name, n_layers, labelpos=None,
                      legpos="lower left", xlabel="$c$",
                      simL=[], save=True):

    fig, axes = plt.subplots(1, 2, sharey=True, figsize=figsize)
    plot_regions(axes[0], name, n_layers, labelpos, legpos, xlabel)

    for L in simL:
        plot_gvc(axes[0], name + "{}_sim".format(L))

    if save:
        fig.tight_layout()
        fig.savefig("Report/multilayer_single_param_{}.pdf".format(name))

    return fig, axes

plot_single_param("ErdosRenyiGraph",
                  [2, 3, 4, 5],
                  labelpos=[(2.25, 0.6), (2.99, 0.55), (3.44, 0.5), (3.65, 0.4)],
                  simL=[1, 2, 3, 4, 5])

plot_single_param("GeometricGraph", [2, 3, 4],
                  labelpos=[(2.25, 0.6), (3.25, 0.5), (4.25, 0.4)],
                  simL=[])

fig, axes = plot_single_param("ScaleFreeGraph", [2],
                          labelpos=[(2.2, 0.6)],
                          xlabel="$\\alpha$",
                          save=False,
                          legpos="lower right",
                          simL=[])

minx = 1.75
ax = axes[0]
ax.set_xlim(minx, ax.get_xlim()[1])

unkown = Rectangle([minx, 0], 2 - minx, 1, color="C1")
ax.add_patch(unkown)
ax.plot([1, 2], [0, 0], linewidth=1, color="C0")
ax.plot([1, 2], [1, 1], linewidth=1, color="C0")

fig.tight_layout()
fig.savefig("Report/multilayer_single_param_ScaleFreeGraph.pdf")
plt.show()
