import json
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

plt.rc("text", usetex=True)
plt.rc("font", size=10, family="palatino linotype")
figsize = (4, 3)

head = "Plot generation/single_param_multiplex/"

def plot_patches(ax, name):
    with open(head + name) as file:
        all_data = json.load(file)

    unkown = []
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
            exist.append(rec)
        elif data["status"] == "unkown":
            unkown.append(rec)

    ax.add_collection(PatchCollection(exist, facecolor="C0"))
    ax.add_collection(PatchCollection(unkown, facecolor="C1"))
    ax.set_xlim(xmin, xmax)

def plot_single_param(name, n_layers, labelpos=None, legpos="upper left"):
    fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize)
    for k, L in enumerate(n_layers):
        fullname = "{}{}.json".format(name, L)
        plot_patches(ax, fullname)
        if labelpos is not None:
            x, y = labelpos[k]
            ax.text(x, y, "$L = {}$".format(L))
                    # bbox=dict(facecolor="white", alpha=0.5, linewidth=0))

    ax.legend((Rectangle((0, 0), 1, 1, facecolor="C0"),
                Rectangle((0, 0), 1, 1, facecolor="C1")),
                ("Solutions exist", "Unkown status"), loc=legpos)
    ax.set_xlabel("$c$")
    ax.set_ylabel("$u$")
    fig.tight_layout()
    fig.savefig("LaTeX/Report/multilayer_single_param_{}.pdf".format(name))

plot_single_param("ErdosRenyiGraph",
                  [2, 3, 4, 5],
                  labelpos=[(2.25, 0.6), (2.99, 0.55), (3.44, 0.5), (3.65, 0.4)])
plt.show()
