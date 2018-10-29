import numpy as np
import matplotlib.pyplot as plt


file = "data.txt"

y = np.loadtxt(file, usecols=0)
plt.plot(y)
plt.show()
