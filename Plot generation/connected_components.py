import networkx as nx
import numpy as np

from matplotlib.colors import Normalize, LogNorm
from matplotlib import cm
from matplotlib import pyplot as plt

c = 1.
n = 100

g = nx.erdos_renyi_graph(n, c/n, seed=2)
comps = [c for c in nx.connected_components(g)]
comps = sorted(comps, key=len)
lens = np.asarray([len(c) for c in comps])
len_set = sorted(list(set(lens)))
cmap = cm.get_cmap("inferno_r")
norm = Normalize(vmin=0, vmax=len(len_set) + 1)
col_dict = {len_set[i]:cmap(norm(i + 1)) for i in range(len(len_set))}

colors = [0 for _ in range(n)]
for (comp, l) in zip(comps, lens):
    for v in comp:
        colors[v] = col_dict[l]

nx.draw(g, pos=nx.spring_layout(g, k=1.5/np.sqrt(n), iterations=100), node_color=colors, node_size=20)
plt.show()
