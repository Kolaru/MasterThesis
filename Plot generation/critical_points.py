import numpy as np
from matplotlib import pyplot as plt

plt.rc("text", usetex=True)
plt.rc("font", size=10, family="palatino linotype")

y_shift = 0.05
umin = 0.6
umax = 1.4
u = np.linspace(umin, umax, 10)
umask = lambda x: np.logical_and(umin < x, x < umax)
capwidth = 0.02

x = np.linspace(0, 2, 10000)
x_low = x[x <= 1]
x_high = x[x > 1]
y = np.zeros_like(x)
y[x > 1] = np.sqrt(x_high - (1 - y_shift)) - np.sqrt(y_shift)
y = 1 - y

ymin = np.min(y) - 0.1*np.abs(np.min(y))
ymax = 1.3

z_low = np.ones_like(x_low)
z_high = 1 - (np.sqrt(x_high - (1 - y_shift)) + 0.05)

fig, axes = plt.subplots(1, 2, figsize=(6.5, 3))

ax = axes[0]

base_curve, = ax.plot(x, np.ones_like(x), color="gray", linestyle="dotted")
ax.plot(x, y, color="gray", linestyle="dotted")
trivial_curve, = ax.plot(x[umask(x)], np.ones_like(x[umask(x)]), linewidth=3)
nontrivial_curve, = ax.plot(x[umask(x)], y[umask(x)], linewidth=3, linestyle="dashed")
ax.plot(u, 1.1*np.ones_like(u), color="black", linewidth=2)
ax.plot([umin, umin], [1.1 - capwidth, 1.1 + capwidth], color="black", linewidth=2)
ax.plot([umax, umax], [1.1 - capwidth, 1.1 + capwidth], color="black", linewidth=2)
ax.plot(1, 1, color="black", marker="o", markersize="8", mec="white", mew="2")

ax.text(1, 1.15, "$U$")

labelpos = (1.7, 0.75)
ax.annotate('Critical point\n $(\lambda^*, u_T)$', xy=(1, 1), xytext=labelpos,
            arrowprops=dict(facecolor='black', shrink=0.15),
            ha="center"
            )

ax.set_ylim(ymin, ymax)
ax.set_xlabel("$\lambda$")
ax.set_ylabel("$u$")
ax.set_xticks([])
ax.set_yticks([])
ax.legend([base_curve, trivial_curve, nontrivial_curve],
          ["All solutions",
           "Trivial implicit function $h_T(\\lambda)$",
           "Non trivial implicit function $h(\\lambda)$"])

ax = axes[1]

base_curve, = ax.plot(x, np.ones_like(x), color="gray", linestyle="dotted")
ax.plot(x_low, z_low, color="gray", linestyle="dotted")
ax.plot(x_high, z_high, color="gray", linestyle="dotted")

trivial_curve, = ax.plot(x[umask(x)], np.ones_like(x[umask(x)]), linewidth=3)
nontrivial_curve, = ax.plot(x_low[umask(x_low)], z_low[umask(x_low)], color="C1", linewidth=3, linestyle="dashed")
ax.plot(x_high[umask(x_high)], z_high[umask(x_high)], color="C1", linewidth=3, linestyle="dashed")
ax.plot(u, 1.1*np.ones_like(u), color="black", linewidth=2)
ax.plot([umin, umin], [1.1 - capwidth, 1.1 + capwidth], color="black", linewidth=2)
ax.plot([umax, umax], [1.1 - capwidth, 1.1 + capwidth], color="black", linewidth=2)
ax.plot(1, 1, color="black", marker="o", markersize="8", mec="white", mew="2")
ax.plot(1, np.max(z_high), color="black", marker="o", markersize="8", mec="white", mew="2")

ax.text(1, 1.15, "$U$")

labelpos = (1.6, 0.6)
ax.annotate('Critical point\n $(\lambda^*, u_T)$', xy=(1, 1), xytext=labelpos,
            arrowprops=dict(facecolor='black', shrink=0.15),
            ha="center"
            )

labelpos = (0.6, 0.4)
ax.annotate('', xy=(1, np.max(z_high)), xytext=labelpos,
            arrowprops=dict(facecolor='black', shrink=0.1),
            ha="center"
            )
ax.text(labelpos[0], labelpos[1] - 0.05, 'Critical point $(\lambda^*, u^\dagger)$',
        ha="center")

ax.set_ylim(ymin, ymax)
ax.set_xlabel("$\lambda$")
ax.set_xticks([])
ax.set_yticks([])
# ax.legend([base_curve, trivial_curve, nontrivial_curve],
#           ["All solutions", "Trivial inverse function $h_T(\lambda)$", "Nontrivial inverse function $h(\lambda)$"],
#           loc="lower left")

fig.tight_layout()
fig.savefig("LaTeX/Report/critical_point.pdf")
plt.show()
