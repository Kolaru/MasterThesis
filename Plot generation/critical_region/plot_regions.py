import json
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

plt.rc("text", usetex=True)
plt.rc("font", size=10, family="palatino linotype")
figsize = (4, 3)

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

    ax.add_collection(PatchCollection(trivial, facecolor="C1", edgecolor="k"))
    ax.add_collection(PatchCollection(nontrivial, facecolor="C0", edgecolor="k"))

def plot_regions(name):
    fullname = "{}.json".format(name,)

    fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize)
    plot_patches(ax, fullname)
    ax.legend((Rectangle((0, 0), 1, 1, facecolor="C0"),
                Rectangle((0, 0), 1, 1, facecolor="C1")),
                ("Non trivial", "Trivial"),
                loc="upper right")

    ax.set_xlabel("$c_1$")
    ax.set_ylabel("$c_2$")
    ax.set_xlim((2, 3))
    ax.set_ylim((2, 3))
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig("Report/critical_region_{}.pdf".format(name))

plot_regions("ErdosRenyiGraph_ErdosRenyiGraph")
plt.show()
