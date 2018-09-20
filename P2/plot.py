import numpy as np
import matplotlib.pyplot as plt

#file = "benchmark.txt"
file = "eigenvec.txt"

x = np.loadtxt(file, usecols=0)
y1 = np.loadtxt(file, usecols=1)
y2 = np.loadtxt(file, usecols=2)

plt.plot(x, y1**2)
plt.plot(x, y2**2)


print(np.trapz(y1**2, x))
#plt.plot(np.log10(n), 3 * np.log10(n) - 5.5)
plt.show()
