import json
import numpy as np

from matplotlib import pyplot as plt

plt.rc("text", usetex=True)
plt.rc("font", size=10, family="palatino linotype")
figsize = (4, 3)

def plot_gcc(name, paramname="c"):
        with open(head + name + ".json") as file:
            data = json.load(file)

        c0 = np.min(data[0]["parameters"])
        c1 = np.max(data[0]["parameters"])

        c = np.linspace(c0, c1, npoints)
        u = np.zeros(npoints)
        for i in range(1000):
            u = g1[name](u, c)

        fig, ax = plt.subplots(figsize=figsize)

        ax.plot(c, 1 - g0[name](u, c), color="k")
        for d, m in zip(data, markers):
            ax.plot(d["parameters"], d["results"], marker=m, linewidth=0,
                label="$n = {}$, $M = {}$".format(d["n"], d["repeat"]))

        ax.set_xlabel("${}$".format(paramname))
        ax.set_ylabel("$S$")
        ax.legend()
        fig.tight_layout()
        fig.savefig("LaTeX/Report/GCC_{}.pdf".format(name))

def geom_g0(z, c):
    return 1/c * z/(1 - (1 - 1/c)*z)

def geom_g1(z, c):
    return 1/c**2 * 1/(1 - (1 - 1/c)*z)**2

def ER_g0(z, c):
    return np.exp(c*(z - 1))

head = "Data/Simulations/"

g0 = {"Geometric": geom_g0, "ER": ER_g0}
g1 = {"Geometric": geom_g1, "ER": ER_g0}

npoints = 1000
markers = "os^d"

plot_gcc("Geometric")
plot_gcc("ER")

plt.show()
