import json
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt

for_ppt = True

if for_ppt:
    from presentation_styling import *
    nodecolor = content
    edgecolor = content
else:
    head = "Report/"
    nodecolor = "black"
    edgecolor = "black"

    def savefig(fig, name):
        fig.savefig(head + name + ".pdf")

RELOAD = True

class RandPos:
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

def norm(vec):
    vec = np.asarray(vec)
    return np.sqrt(np.sum(vec**2))

def load_network(name):
    return nx.read_edgelist("Data/real-networks/{}/out.{}".format(name, name),
        comments="%")

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

def swap(d, k1, k2):
    k1 = str(k1)
    k2 = str(k2)
    temp = d[k1]
    d[k1] = d[k2]
    d[k2] = temp

def is_valid(center, excluded_pos, tol, radius):
    if norm(center) < radius:
        return False

    for pos in excluded_pos:
        if norm(pos - center) < tol:
            return False

    return True

def network_subplots():
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.axis("off")
    ax.set_aspect("equal")
    return fig, ax

nodesize = 2

radius = 0.4
np.random.seed(0)
randpos = RandPos(radius + 0.02, 1)

name = "subelj_euroroad"
g = load_network(name)
comps = sorted(list(nx.connected_component_subgraphs(g)), key=len, reverse=True)
gcc = comps[0]
pos = node_positions(g, name, nx.kamada_kawai_layout)

fig, ax = network_subplots()
nx.draw_networkx_edges(gcc, pos=pos, edge_color=edgecolor, alpha=0.2)
nx.draw_networkx_nodes(gcc, pos=pos, node_size=nodesize, node_color=nodecolor)

excluded_pos = [pos[k] for k in gcc if norm(pos[k]) > radius]
for r, comp in enumerate(comps[1:]):
    valid_pos = False
    ntest = 0
    while not valid_pos:
        center = randpos.next()
        if is_valid(center, excluded_pos, 0.05*np.sqrt(len(comp)), radius):
            excluded_pos.append(center)
            valid_pos = True

    comp_pos = nx.kamada_kawai_layout(comp, center=center, scale=0.03*np.sqrt(len(comp)))
    nx.draw_networkx_edges(comp, pos=comp_pos, edge_color=edgecolor, alpha=0.2)
    nx.draw_networkx_nodes(comp, pos=comp_pos, node_size=nodesize, node_color=nodecolor)

fig.tight_layout()
savefig(fig, "network-{}".format(name))


name = "arenas-jazz"
g = load_network(name)
pos = node_positions(g, name, nx.kamada_kawai_layout)

pos["25"] = np.asarray([-0.55, -0.1])
pos["195"] = np.asarray([0.55, 0.35])
pos["175"] = np.asarray([0.5, 0.3])
pos["162"] = np.asarray([0.5, 0.])
pos["79"] = np.asarray([0.4, -0.7])
pos["78"] = np.asarray([0.5, -0.6])
pos["77"] = np.asarray([0.5, -0.75])
pos["198"] = np.asarray([-0.5, -0.3])

fig, ax = network_subplots()
# nx.draw_networkx_labels(g, pos=pos)
nx.draw_networkx_edges(g, ax=ax, pos=pos, edge_color=edgecolor, alpha=0.2)
nx.draw_networkx_nodes(g, ax=ax, pos=pos, node_size=nodesize, node_color=nodecolor)
fig.tight_layout()
savefig(fig, "network-{}".format(name))

name = "US-power-grid"
g = load_network(name)
pos = node_positions(g, name, nx.kamada_kawai_layout)

fig, ax = network_subplots()
nx.draw_networkx_edges(g, ax=ax, pos=pos, edge_color=edgecolor, alpha=0.2)
nx.draw_networkx_nodes(g, ax=ax, pos=pos, node_size=nodesize, node_color=nodecolor)
fig.tight_layout()
savefig(fig, "network-{}".format(name))


radius = 0.75
np.random.seed(2)
randpos = RandPos(radius + 0.02, 1.1)

name = "maayan-vidal"
g = load_network(name)
comps = sorted(list(nx.connected_component_subgraphs(g)), key=len, reverse=True)
gcc = comps[0]
pos = node_positions(g, name, nx.kamada_kawai_layout, generate_new=False)

fig, ax = network_subplots()
nx.draw_networkx_edges(gcc, pos=pos, edge_color=edgecolor, alpha=0.2)
nx.draw_networkx_nodes(gcc, pos=pos, node_size=nodesize, node_color=nodecolor)

excluded_pos = [pos[k] for k in gcc if norm(pos[k]) > radius]
for r, comp in enumerate(comps[1:]):
    valid_pos = False
    ntest = 0
    while not valid_pos:
        center = randpos.next()
        if is_valid(center, excluded_pos, 0.05*len(comp), radius):
            excluded_pos.append(center)
            valid_pos = True

    comp_pos = nx.spring_layout(comp, random_state=r, center=center, scale=0.02)
    nx.draw_networkx_edges(comp, pos=comp_pos, edge_color=edgecolor, alpha=0.2)
    nx.draw_networkx_nodes(comp, pos=comp_pos, node_size=nodesize, node_color=nodecolor)

fig.tight_layout()
savefig(fig, "network-{}".format(name))
plt.show()
