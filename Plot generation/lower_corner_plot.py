import json
import numpy as np
import matplotlib

matplotlib.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

plt.rc("text", usetex=True)
plt.rc("font", size=10, family="palatino linotype")
figsize = (4.5, 3.5)

capsize = 0.04

xmin = -0.75
xmax = -0.5
xx = np.linspace(xmin, xmax, 1000)

fig, ax = plt.subplots(figsize=figsize)

# Dummy function
def f(x):
    return x**2

# Solution curve
solnt, = ax.plot(xx, f(xx), color="C3")
soltriv, = ax.plot([xmin, xmax], [1, 1], color="C0")

# lambda line
capsize = 0.04
lambdalo = (xmax - xmin)/3 + xmin
lambdahi = 2*(xmax - xmin)/3 + xmin
lambdapos = 1.07
ax.plot([lambdalo, lambdahi], [lambdapos, lambdapos], color="k", lw=2)
ax.plot([lambdalo, lambdalo], [lambdapos - capsize/2, lambdapos + capsize/2], color="k", lw=2)
ax.plot([lambdahi, lambdahi], [lambdapos - capsize/2, lambdapos + capsize/2], color="k", lw=2)
ax.text( (lambdalo + lambdahi)/2, lambdapos + 0.04, "$\\boldsymbol{\\Lambda}$", ha="center")

# Z box
capsize = 0.008
zminx = lambdalo
zmaxx = lambdahi
zminy = f(lambdahi) - 0.02
zmaxy = 1

zw = zmaxx - zminx
zh = zmaxy - zminy

recz = Rectangle([zminx, zminy], zw, zh, color="C1", alpha=0.2)
ax.add_patch(recz)

# Z line
zpos = lambdahi + 0.011
ax.plot([zpos, zpos], [zminy, zmaxy], color="k", lw=2)
ax.plot([zpos - capsize/2, zpos + capsize/2], [zminy, zminy], color="k", lw=2)
ax.plot([zpos - capsize/2, zpos + capsize/2], [zmaxy, zmaxy], color="k", lw=2)
ax.text(zpos + 0.008, (zminy + zmaxy)/2, "$\\boldsymbol{\\mathrm{Z}}$", va="center")

# Corner line
cpos = lambdalo - 0.011
clo = zminy
chi = 0.25 + 0.75*zminy
ax.plot([cpos, cpos], [clo, chi], color="k", lw=2)
ax.plot([cpos - capsize/2, cpos + capsize/2], [clo, clo], color="k", lw=2)
ax.plot([cpos - capsize/2, cpos + capsize/2], [chi, chi], color="k", lw=2)
ax.text(cpos - 0.03, (clo + chi)/2, "$\\boldsymbol{\\mathrm{C}}^{(\\gamma)}$", va="center")

reccorner = Rectangle([zminx, clo], zw, (chi - clo), lw=2.5, edgecolor="C1", facecolor="none", alpha=1)
ax.add_patch(reccorner)


ax.set_xticks([])
ax.set_yticks([])
ax.set_ylim(-0.1, 1.2)

ax.set_xlabel("$\\lambda$")
ax.set_ylabel("$z$")

ax.legend((soltriv, solnt, recz, reccorner),
          ["Trivial solutions", "Non trivial solutions",
          "Refined region", "Lower $\\gamma$ corner"],
          ncol=2)
fig.savefig("Report/lower_corner.pdf")
plt.show()
