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

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize)
plot_patches(ax, "ErdosRenyiGraph2.json")
plot_patches(ax, "ErdosRenyiGraph3.json")
plot_patches(ax, "ErdosRenyiGraph4.json")
plot_patches(ax, "ErdosRenyiGraph5.json")

ax.legend((Rectangle((0, 0), 1, 1, facecolor="C0"),
            Rectangle((0, 0), 1, 1, facecolor="C1")),
            ("Solutions exist", "Unkown status"), loc="upper left")
ax.set_xlabel("$c$")
ax.set_ylabel("$u$")
fig.tight_layout()
plt.show()
