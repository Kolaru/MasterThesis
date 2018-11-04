import json
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

plt.rc("text", usetex=True)
plt.rc("font", size=10, family="palatino linotype")
figsize = (3, 3)

head = "Plot generation/critical_region/"

def to_rect(data):
    recs = []
    for bounds in data:
        yb, xb = bounds
        w = xb[1] - xb[0]
        h = yb[1] - yb[0]
        recs.append(Rectangle([xb[0], yb[0]], w, h))
    return recs

def plot_patches(ax, name):
    with open(head + name) as file:
        all_data = json.load(file)

    trivial = to_rect(all_data["trivial"])
    nontrivial = to_rect(all_data["nontrivial"])

    ax.add_collection(PatchCollection(trivial, facecolor="C1"))
    ax.add_collection(PatchCollection(trivial, facecolor="none", edgecolor="white", alpha=0.2))
    ax.add_collection(PatchCollection(nontrivial, facecolor="C0"))
    ax.add_collection(PatchCollection(nontrivial, facecolor="none", edgecolor="white", alpha=0.2))

def plot_regions(name, legpos="upper right",
                 insetxlim=(1, 2), insetylim=(1, 2),
                 insetlocs=(2, 4)):
    fullname = "{}.json".format(name)
    numname = "{}_numerical.json".format(name)

    fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize)

    plot_patches(ax, fullname)

    with open(head + numname) as file:
        numerical_data = json.load(file)

    ax.plot(numerical_data["x"], numerical_data["y"], color="k")

    if legpos is not None:
        ax.legend((Rectangle((0, 0), 1, 1, facecolor="C0"),
                    Rectangle((0, 0), 1, 1, facecolor="C1")),
                    ("Non trivial", "Trivial"),
                    loc=legpos)


    ax.set_xlabel("$c_1$")
    ax.set_ylabel("$c_2$")
    ax.set_xlim((2, 3))
    ax.set_ylim((2, 3))
    ax.set_aspect("equal")

    axins = inset_axes(ax, width="35%", height="35%", loc="lower left")
    plot_patches(axins, fullname)
    axins.plot(numerical_data["x"], numerical_data["y"], color="k", lw=2)

    axins.set_xlim(insetxlim)
    axins.set_ylim(insetylim)
    mark_inset(ax, axins, loc1=insetlocs[0], loc2=insetlocs[1], ec="k", linewidth="1.")
    axins.set_xticks([])
    axins.set_yticks([])

    fig.tight_layout()
    fig.savefig("Report/critical_region_{}.pdf".format(name))

plot_regions("ErdosRenyiGraph_ErdosRenyiGraph",
             insetxlim=(2.15, 2.25),
             insetylim=(2.74, 2.84),
             insetlocs=(1, 2))
plot_regions("GeometricGraph_GeometricGraph", legpos=None,
             insetxlim=(2.39, 2.49),
             insetylim=(2.76, 2.86))
plt.show()
