from matplotlib import pyplot as plt

import numpy as np

plt.rc("text", usetex=True)

def g1(z):
    return np.exp(1.3*(z - 1))

zz = np.linspace(0, 1, 10000)

fig, ax = plt.subplots()

zk = [0.]
xk = []
yk = []

for _ in range(10):
    a = zk[-1]
    xk.extend([a, a])
    yk.append(a)
    zk.append(g1(a))
    yk.append(zk[-1])

u = 0

for _ in range(1000):
    u = g1(u)

g1line, = ax.plot(zz, g1(zz), label="$g_1(z)$")
idline, = ax.plot(zz, zz, label="$z = y$")
iterline, = ax.plot(xk, yk, "-", color="gray", label="Fixpoint iteration")
ax.plot(u, u, "o", color="gray")
ax.annotate("$u = g_1(u)$", xy=(u, u), xytext=(u + 0.2, u - 0.2),
                arrowprops=dict(facecolor='black', shrink=0.2),
                ha="right"
                )

ax.set_xlabel("$z$")
ax.set_ylabel("$y$")
ax.set_aspect("equal")
ax.legend()
fig.tight_layout()
fig.savefig("LaTeX/Report/u_solution_graphically.pdf")
plt.show()
