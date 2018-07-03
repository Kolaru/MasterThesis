import numpy as np
from matplotlib import pyplot as plt

plt.rc("text", usetex=True)

continuous = False
discontinuous = True

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

z_low = np.ones_like(x_low)
z_high = 1 - (np.sqrt(x_high - (1 - y_shift)) + 0.05)

if continuous:
    fig, ax = plt.subplots()
    base_curve, = ax.plot(x, np.ones_like(x), color="gray", linestyle="dotted")
    ax.plot(x, y, color="gray", linestyle="dotted")
    trivial_curve, = ax.plot(x[umask(x)], np.ones_like(x[umask(x)]), linewidth=3)
    nontrivial_curve, = ax.plot(x[umask(x)], y[umask(x)], linewidth=3, linestyle="dashed")
    ax.plot(u, 1.1*np.ones_like(u), color="black", linewidth=2)
    ax.plot([umin, umin], [1.1 - capwidth, 1.1 + capwidth], color="black", linewidth=2)
    ax.plot([umax, umax], [1.1 - capwidth, 1.1 + capwidth], color="black", linewidth=2)
    ax.plot(1, 1, color="black", marker="o", markersize="8", mec="white", mew="2")

    ax.text(1, 1.15, "$U$")
    ax.annotate('Critical point $(\lambda^*, u_T)$', xy=(1, 1), xytext=(0.8, 0.8),
                arrowprops=dict(facecolor='black', shrink=0.1),
                ha="right"
                )

    ax.set_ylim(np.min(y) - 0.1*np.abs(np.min(y)), 1.5)
    ax.set_xlabel("$\lambda$")
    ax.set_ylabel("$u$")
    ax.legend([base_curve, trivial_curve, nontrivial_curve],
              ["All solutions", "Trivial inverse function $h_T(\lambda)$", "Nontrivial inverse function $h(\lambda)$"],
              loc="lower left")

if discontinuous:
    fig, ax = plt.subplots()

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
    ax.annotate('Critical point $(\lambda^*, u_T)$', xy=(1, 1), xytext=(0.8, 0.8),
                arrowprops=dict(facecolor='black', shrink=0.1),
                ha="right"
                )

    ax.annotate('Critical point $(\lambda^*, u^\dagger)$', xy=(1, np.max(z_high)), xytext=(1.3, 0.8),
                arrowprops=dict(facecolor='black', shrink=0.1),
                ha="left"
                )

    ax.set_ylim(np.min(y) - 0.1*np.abs(np.min(y)), 1.5)
    ax.set_xlabel("$\lambda$")
    ax.set_ylabel("$u$")
    ax.legend([base_curve, trivial_curve, nontrivial_curve],
              ["All solutions", "Trivial inverse function $h_T(\lambda)$", "Nontrivial inverse function $h(\lambda)$"],
              loc="lower left")

plt.show()
