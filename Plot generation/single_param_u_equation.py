import json
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

plt.rc("text", usetex=True)
plt.rc("font", size=10, family="palatino linotype")
figsize = (4, 3)

def g0(z):
    return np.exp(2.6*(z - 1))

def ufunc(z):
    return 1 - (1 - g0(z))**2

fig, ax = plt.subplots(figsize=(3.5, 3.5))

zz = np.linspace(0, 1, 10000)
g1line, = ax.plot(zz, ufunc(zz), label="$\\psi(z)$")
idline, = ax.plot(zz, zz, label="$z = y$")
# ax.plot(u, u, "o", color="gray")
# ax.annotate("$u = g_1(u)$", xy=(u, u), xytext=(u + 0.3, u - 0.3),
#                 arrowprops=dict(facecolor='black', shrink=0.2),
#                 ha="right"
#                 )

ax.set_xlabel("$z$")
ax.set_ylabel("$y$")
ax.set_aspect("equal")
ax.legend()
fig.tight_layout()
fig.savefig("LaTeX/Report/single_param_u_solution_graphically.pdf")
plt.show()
