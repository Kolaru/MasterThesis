import json
import mpmath as mp
import numpy as np

from matplotlib import pyplot as plt
from scipy.special import zeta

polylog = np.frompyfunc(lambda a, z: np.real(mp.fp.polylog(a, z)), 2, 1)

plt.rc("text", usetex=True)
plt.rc("font", size=10, family="palatino linotype")
figsize = (6.6, 2.6)

for_ppt = True

if for_ppt:
    import sys
    sys.path.append("Plot generation")

    from presentation_styling import *
    bgcolor = bg
    linecolor = content
else:
    head = "Report/"
    bgcolor = "white"
    linecolor = "k"

    def savefig(fig, name):
        fig.savefig(head + name + ".pdf")


def plot_gcc(name, paramname="c", critpoint=1, critlabely=0.3):
        with open(head + name + ".json") as file:
            data = json.load(file)

        c0 = np.min(data[0]["parameters"])
        c1 = np.max(data[0]["parameters"])

        c = np.linspace(c0 - 0.15, c1 + 0.15, npoints)
        u = np.zeros(npoints)
        for i in range(1000):
            u = g1[name](u, c)

        fig, axes = plt.subplots(1, 3, figsize=figsize, sharey=True)

        ax = axes[0]

        for ax in axes:
            ax.axvline(critpoint, color="gray", zorder=-100)
            ax.text(critpoint, critlabely, "$\\mathcal{R}$", ha="center",
                    bbox=dict(facecolor=bgcolor, linewidth=0))

            ax.set_xlabel("${}$".format(paramname))
            ax.set_xlim((c0 - 0.1, c1 + 0.1))

            ax.plot(c, 1 - g0[name](u, c), color=linecolor, zorder=-5)

        for k, (d, m) in enumerate(zip(data, markers)):
            axes[k].errorbar(d["parameters"], d["sizes"], yerr=d["stds"],
                             marker=".",
                             color="C{}".format(k),
                             linewidth=0,
                             capsize=4,
                             elinewidth=2,
                             markeredgewidth=2)
            axes[k].set_title(label="$n = {}$, $M = {}$".format(d["n"], d["repeat"]))

        axes[0].set_ylabel("$S$")
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.2, top=0.85, wspace=0.1)
        savefig(fig, "GCC_{}".format(name))

def geom_g0(z, c):
    return 1/c * z/(1 - (1 - 1/c)*z)

def geom_g1(z, c):
    return 1/c**2 * 1/(1 - (1 - 1/c)*z)**2

def ER_g0(z, c):
    return np.exp(c*(z - 1))

def sf_g0(z, a):
    return polylog(a, z)/zeta(a)

def sf_g1(z, a):
    res = np.zeros_like(z)
    res[z != 0] = polylog(a[z != 0]-1, z[z != 0])/(z[z != 0] * zeta(a[z != 0]-1))
    res[z == 0] = 1/(2**a[z == 0])
    return res

head = "Plot generation/gcc_plots/"

g0 = {"Geometric": geom_g0, "ER": ER_g0, "Scalefree": sf_g0}
g1 = {"Geometric": geom_g1, "ER": ER_g0, "Scalefree": sf_g1}

npoints = 1000
markers = "os^d"

plot_gcc("Geometric", critpoint=1.5, critlabely=0.7)
plot_gcc("ER", critpoint=1, critlabely=0.5)
plot_gcc("Scalefree", critpoint=3.47875, paramname="\\alpha", critlabely=0.5)

plt.show()
