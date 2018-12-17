import json
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from presentation_styling import *

nodecolor = content
edgecolor = content
nodesize = 2


class DiskRand:
    def __init__(self, radmin, radmax, N=1000):
        self.radmin = radmin
        self.radmax = radmax
        self.N = N
        self.new_precalc()

    def new_precalc(self):
        self.k = 0
        randr = np.random.uniform(self.radmin, self.radmax, size=self.N)
        randphi = np.random.uniform(0, 2*np.pi, size=self.N)
        self.randpos = np.zeros((self.N, 2))
        self.randpos[:, 0] = randr * np.cos(randphi)
        self.randpos[:, 1] = randr * np.sin(randphi)

    def next(self):
        if self.k == self.N:
            self.new_precalc()

        res = self.randpos[self.k]
        self.k += 1
        return res


def node_positions(g, name, layout, generate_new=False):
    if generate_new:
        pos = layout(g)
        for key in pos:
            pos[key] = pos[key].tolist()
        with open("Plot generation/layouts/{}.json".format(name), "w") as file:
            json.dump(pos, file)
    else:
        with open("Plot generation/layouts/{}.json".format(name), "r") as file:
            pos = json.load(file)

    return pos


def norm(vec):
    vec = np.asarray(vec)
    return np.sqrt(np.sum(vec**2))


def is_valid(center, excluded_pos, tol, radius):
    if norm(center) < radius:
        return False

    for pos in excluded_pos:
        if norm(pos - center) < tol:
            return False

    return True


np.random.seed(2)

def find_nearest3(i, allpos):
    v = allpos[i]
    min1 = 10
    min2 = 10
    min3 = 10
    k1 = -1
    k2 = -1
    k3 = -1
    for k, w in enumerate(allpos):
        if k != i:
            nw = norm(w - v)
            if nw < min1:
                min3 = min2
                k3 = k2
                min2 = min1
                k2 = k1
                min1 = nw
                k1 = k
            elif nw < min2:
                min3 = min2
                k3 = k2
                min2 = nw
                k2 = k
            elif nw < min3:
                min3 = nw
                k3 = k

    return k1, k2, k3


def in_gcc(phi):
    return np.pi/3 < phi <= np.pi


def update(frac, ringw, R, rho, allpos):
    r = frac * R + ringw
    dA = np.pi*(r**2 - (r - ringw)**2)
    N = int(rho * dA)
    nprev = len(allpos)

    if N > 0:
        ax.set_xlim(-r, r)
        ax.set_ylim(-r, r)

        L = np.random.uniform(r - ringw, r, N)
        phi = np.random.uniform(0, 2*np.pi, N)

        pos = L*np.asarray([np.cos(phi), np.sin(phi)])
        allpos.extend(list(pos.T))
        ax.plot(*pos, color=content, marker=".", markersize=1, linewidth=0)
        edges = []

        pos = pos.T

        for k in range(len(pos)):
            v = pos[k]
            k1, k2, k3 = find_nearest3(nprev + k, allpos)

            if in_gcc(phi[k]):
                if k1 > -1:
                    w = allpos[k1]
                    ax.plot([v[0], w[0]], [v[1], w[1]], linewidth=1, color=content, alpha=alpha)

                if k2 > -1:
                    w = allpos[k2]
                    ax.plot([v[0], w[0]], [v[1], w[1]], linewidth=1, color=content, alpha=alpha)

                if k3 > -1:
                    w = allpos[k3]
                    ax.plot([v[0], w[0]], [v[1], w[1]], linewidth=1, color=content, alpha=alpha)
            else:
                if np.random.rand() > 0.7:
                    w = allpos[k1]
                    ax.plot([v[0], w[0]], [v[1], w[1]], linewidth=1, color=content, alpha=alpha)

fig, ax = plt.subplots()
ax.set_aspect("equal")
ax.axis("off")

alpha = 0.5

R = 1
T = 5
interval = 30  # ms
df = interval/1000 / T
fracs = np.arange(0, 1, df)
rho = 200/(np.pi*R**2)
minr = 0.05
ringw = 0.05
allpos = []

anim = FuncAnimation(fig, update, frames=fracs, repeat=False, interval=interval,
                     fargs=(ringw, R, rho, allpos))

# nx.draw_networkx_edges(g, ax=ax, pos=pos, edge_color=edgecolor, alpha=alpha)
# nx.draw_networkx_nodes(g, ax=ax, pos=pos, node_size=nodesize, node_color=nodecolor)

anim.save("Presentation/network_growth.avi", savefig_kwargs=dict(facecolor=bg))

plt.show()
