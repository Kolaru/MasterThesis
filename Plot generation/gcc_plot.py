from pandas import read_json
from matplotlib import pyplot as plt
import numpy as np

data = read_json("Data/Simulations/Erdos_Renyi.json")
plt.plot(data["parameters"][0], data["results"][0], "o")
plt.show()
