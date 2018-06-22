from matplotlib import pyplot as plt
import numpy as np

plt.rc("text", usetex=True)

# Poisson degree dist with p_0 = 0
class PoissonClass:
    def pks(self, c, n):
        pks = np.zeros(n)
        pk = 1
        for k in np.arange(1, n+1):
            pk *= c/(k)
            pks[k-1] = pk

        return np.exp(-c)*pks/(1 - np.exp(-c))

    def g0(self, c, z):
        return np.exp(c*(z-1)) - np.exp(-c)

    def g1(self, c, z):
        return np.exp(c*(z-1))

Poisson = PoissonClass()

n = 25
cc = [1.1, 1.25, 1.7]
markers = ["o", "s", "d", "^"]
for c, m in zip(cc, markers):
    u = 0
    for _ in range(1000):
        u = Poisson.g1(c, u)

    S = 1 - Poisson.g0(c, u)
    drifts = 1 - u**np.arange(1, n+1)
    plt.plot(drifts/S, m)

plt.legend(["$c = {}$".format(c) for c in cc])
plt.xlabel("$k$")
plt.ylabel("$(1 - u^k)/S$")
plt.show()
