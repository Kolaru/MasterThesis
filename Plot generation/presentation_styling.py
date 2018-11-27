from matplotlib import pyplot as plt

head = "Presentation/"

def savefig(fig, name):
    fig.savefig(head + name + ".png", transparent=True)

bg = "#002b36"
bghl = "#073642"

lowest = "#586e75"
lower = "#657b83"
low = "#839496"

secondary = "#93a1a1"
content = "#eee8d5"
emphasized = "#fdf6e3"

plt.rc("axes", facecolor=bg,
               edgecolor=content)


plt.rc("figure", facecolor=bg,
               edgecolor=bg)


plt.rc("savefig", facecolor=bg,
               edgecolor=bg)
