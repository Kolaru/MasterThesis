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
                 insetlocs=(2, 4), lims=(2, 3),
                 insetpos="lower left",
                 xlabel="$c_1$", ylabel="$c_2$",
                 numerical=True):
    fullname = "{}.json".format(name)
    numname = "{}_numerical.json".format(name)

    fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize)

    plot_patches(ax, fullname)

    if legpos is not None:
        ax.legend((Rectangle((0, 0), 1, 1, facecolor="C0"),
                    Rectangle((0, 0), 1, 1, facecolor="C1")),
                    ("Non trivial", "Trivial"),
                    loc=legpos)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect("equal")

    axins = inset_axes(ax, width="35%", height="35%", loc=insetpos)
    plot_patches(axins, fullname)

    axins.set_xlim(insetxlim)
    axins.set_ylim(insetylim)
    mark_inset(ax, axins, loc1=insetlocs[0], loc2=insetlocs[1], ec="k", linewidth="1.")
    axins.set_xticks([])
    axins.set_yticks([])

    if numerical:
        with open(head + numname) as file:
            numerical_data = json.load(file)

        ax.plot(numerical_data["x"], numerical_data["y"], color="k")
        axins.plot(numerical_data["x"], numerical_data["y"], color="k", lw=2)

    fig.tight_layout()
    fig.savefig("Report/critical_region_{}.pdf".format(name))

plot_regions("ErdosRenyiGraph_ErdosRenyiGraph",
             insetxlim=(2.15, 2.25),
             insetylim=(2.74, 2.84),
             insetlocs=(1, 2))
plot_regions("GeometricGraph_GeometricGraph", legpos=None,
             insetxlim=(2.39, 2.49),
             insetylim=(2.76, 2.86))

plot_regions("ScaleFreeGraph_ScaleFreeGraph", legpos=None,
             insetxlim=(2.14, 2.18),
             insetylim=(2.27, 2.31),
             insetpos="upper right",
             lims=(2.1, 2.5),
             xlabel="$\\alpha_1$", ylabel="$\\alpha_2$",
             numerical=True)
plt.show()
