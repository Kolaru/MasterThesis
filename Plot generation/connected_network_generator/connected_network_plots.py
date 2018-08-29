import json
import numpy as np

from matplotlib import pyplot as plt

plt.rc("text", usetex=True)
plt.rc("font", size=14)

head = "Plot generation/connected_network_generator/"

with open(head + "ER_verif.json") as file:
    all_data = json.load(file)

fig, ax = plt.subplots()
K = 20
ks = np.arange(K)
smooth_ks = np.linspace(0, K - 1, 1000)

plots = []
markers = "o^sdxd"

for i, data in enumerate(all_data):
    col = "C" + str(i)
    m = markers[i]
    c = data["c"]
    u = data["u"]
    urecon = data["urecon"]
    pk = np.asarray(data["pk"])
    rk = np.asarray(data["rk"])

    print(u - urecon)

    S = 1 - (np.exp(c*(u - 1)) - np.exp(-c))/(1 - np.exp(-c))

    p, = ax.plot(ks + 1, rk[ks]/pk[ks], m, color=col, label="$c={}$".format(c), markersize=8)
    ax.plot(smooth_ks + 1, (1 - u**(smooth_ks+1))/S, color=col)

    plots.append(p)

ax.set_xlabel("$k$")
ax.set_ylabel("$r_k/p_k$")
ax.set_xticks([1, 5, 10, 15, 20])
ax.legend(handles=plots)

plt.show()
