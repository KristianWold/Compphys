import numpy as np
import matplotlib.pyplot as plt

#file = "benchmark.txt"
file = "eigenvec.txt"

x = np.loadtxt(file, usecols=0)
y1 = np.loadtxt(file, usecols=1)
y2 = np.loadtxt(file, usecols=2)
y3 = np.loadtxt(file, usecols=3)
y4 = np.loadtxt(file, usecols=4)

plt.plot(x, y1**2)
plt.plot(x, y2**2)
plt.plot(x, y3**2)
plt.plot(x, y4**2)
plt.xlabel("|x1 - x2|")
plt.ylabel("probability")
plt.legend(["non-interacting", "Coulomb repulsion, w = 1",
            "Coulomb repulsion, w = 0.3", "Coulomb repulsion, w = 0.06"])

print(np.trapz(y1**2, x))
#plt.plot(np.log10(n), 3 * np.log10(n) - 5.5)
plt.show()
