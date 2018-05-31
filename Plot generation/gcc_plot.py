from pandas import read_json
from matplotlib import pyplot as plt

data = read_json("Data/Erdos_Renyi.json")
plt.plot(data["parameters"][0], data["results"][0])
plt.show()
