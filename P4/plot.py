import numpy as np
import matplotlib.pyplot as plt

file = "results/data.dat"


E = np.fromfile(file, dtype="int32", count=-1)

plt.plot(E)
plt.show()
