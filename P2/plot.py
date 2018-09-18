import numpy as np
import matplotlib.pyplot as plt

file = "benchmark.txt"

n = np.loadtxt(file, usecols=0)
t = np.loadtxt(file, usecols=1)

plt.plot(np.log10(n), np.log10(t))
plt.plot(np.log10(n), 3 * np.log10(n) - 5.5)
plt.show()
