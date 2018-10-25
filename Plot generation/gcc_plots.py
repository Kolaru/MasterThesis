import json
import mpmath as mp
import numpy as np

from matplotlib import pyplot as plt
from scipy.special import zeta

polylog = np.frompyfunc(lambda a, z: np.real(mp.fp.polylog(a, z)), 2, 1)

plt.rc("text", usetex=True)
plt.rc("font", size=10, family="palatino linotype")
figsize = (4, 3)

def plot_gcc(name, paramname="c", critpoint=1, critlabely=0.3):
        with open(head + name + ".json") as file:
            data = json.load(file)

        c0 = np.min(data[0]["parameters"])
        c1 = np.max(data[0]["parameters"])

        c = np.linspace(c0, c1, npoints)
        u = np.zeros(npoints)
        for i in range(1000):
            u = g1[name](u, c)

        fig, ax = plt.subplots(figsize=figsize)

        ax.axvline(critpoint, color="gray")
        ax.text(critpoint, critlabely, "$\\mathcal{R}$", ha="center",
                bbox=dict(facecolor="white", linewidth=0))

        ax.plot(c, 1 - g0[name](u, c), color="k")
        for d, m in zip(data, markers):
            ax.plot(d["parameters"], d["results"], marker=m, linewidth=0,
                label="$n = {}$, $M = {}$".format(d["n"], d["repeat"]))

        ax.set_xlabel("${}$".format(paramname))
        ax.set_ylabel("$S$")
        ax.legend()
        fig.tight_layout()
        fig.savefig("Report/GCC_{}.pdf".format(name))

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

print(sf_g1(np.asarray([0]), np.asarray([2.5])))

head = "Data/Simulations/"

g0 = {"Geometric": geom_g0, "ER": ER_g0, "Scalefree": sf_g0}
g1 = {"Geometric": geom_g1, "ER": ER_g0, "Scalefree": sf_g1}

npoints = 1000
markers = "os^d"

plot_gcc("Geometric", critpoint=1.5, critlabely=0.4)
plot_gcc("ER", critpoint=1)
plot_gcc("Scalefree", critpoint=3.47875, paramname="\\alpha")

plt.show()
