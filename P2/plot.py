import numpy as np
import matplotlib.pyplot as plt
"""
file = "benchmark.txt"

n = np.log10(np.loadtxt(file, usecols=0))
t1 = np.log10(np.loadtxt(file, usecols=1))
t2 = np.log10(np.loadtxt(file, usecols=2))

fit1 = np.poly1d(np.polyfit(n, t1, 1))
fit2 = np.poly1d(np.polyfit(n, t2, 1))

plt.plot(n, t1)
plt.plot(n, fit1(n), ":")
plt.plot(n, t2)
plt.plot(n, fit2(n), ":")
plt.xlabel("log10(n)")
plt.xlabel("log10(t)")
plt.legend(["Jacobi", fit1, "Arma", fit2])
plt.show()
"""

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

plt.show()
