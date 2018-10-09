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

x1 = np.loadtxt(file, usecols=0)
y1 = np.loadtxt(file, usecols=1)

#x2 = np.loadtxt(file, usecols=6)
#y2 = np.loadtxt(file, usecols=7)

#x3 = np.loadtxt(file, usecols=12)
#y3 = np.loadtxt(file, usecols=13)

plt.gca().set_aspect("equal")

plt.plot(x1, y1)
#plt.plot(x2, y2)
#plt.plot(x3, y3)
plt.show()
