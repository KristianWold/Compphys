import numpy as np
import matplotlib.pyplot as plt

# Set fontsizes in figures
params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large'}
plt.rcParams.update(params)

file = "data.txt"

x = np.loadtxt(file, usecols=0)
y = np.loadtxt(file, usecols=1)

plt.gca().set_aspect("equal")

plt.plot(x, y)
plt.show()
