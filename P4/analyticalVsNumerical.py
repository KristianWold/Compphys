import numpy as np
import matplotlib.pyplot as plt

file = "results/evolution_L=2.txt"
t = np.loadtxt(file, usecols=0)
T = np.linspace(t[0], t[-1], 1000)
E_anal = -(8 * np.sinh(8 / T)) / (np.cosh(8 / T) + 3) / 4
M_anal = (2 * np.exp(8 / T) + 4) / (np.cosh(8 / T) + 3) / 4
Cv_anal = 1 / T**2 * 64 * (3 * np.cosh(8 / T) + 1) / \
    (np.cosh(8 / T) + 3)**2 / 4
X_anal = 1 / T * ((8 * np.exp(8 / T) + 8) /
                  (np.cosh(8 / T) + 3) - (4 * M_anal)**2) / 4


E = np.loadtxt(file, usecols=1)
M = np.loadtxt(file, usecols=2)
Cv = np.loadtxt(file, usecols=3)
X = np.loadtxt(file, usecols=4)

fig = plt.figure()
plt.subplot(2, 2, 1)
plt.xlabel("T")
plt.ylabel("<E>")
plt.plot(t, E, linewidth=0.5)
plt.plot(T, E_anal, linewidth=0.5)
plt.grid()

plt.subplot(2, 2, 2)
plt.xlabel("T")
plt.ylabel("<|M|>")
plt.plot(t, M, linewidth=0.5)
plt.plot(T, M_anal, linewidth=0.5)
plt.grid()

plt.subplot(2, 2, 3)
plt.xlabel("T")
plt.ylabel("<Cv>")
plt.plot(t, Cv, linewidth=0.5)
plt.plot(T, Cv_anal, linewidth=0.5)
plt.grid()

plt.subplot(2, 2, 4)
plt.xlabel("T")
plt.ylabel("X")
plt.plot(t, X, linewidth=0.5)
plt.plot(T, X_anal, linewidth=0.5)
plt.grid()

fig.savefig("plots/numericalVsAnalytical.pdf")
