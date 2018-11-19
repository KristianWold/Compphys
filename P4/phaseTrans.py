import os
import numpy as np
import matplotlib.pyplot as plt
import math as m
import sys
import scipy.stats as stats

# phase transtion
# -----------------------------------------------------------------------------
L = np.array([40, 60, 100, 160])
t = np.loadtxt("results/evolution_L=40.txt", usecols=0)

E = []
M = []
Cv = []
X = []

for l in L:
    file = "results/evolution_L=%s.txt" % (l)
    E.append(np.loadtxt(file, usecols=1))
    M.append(np.loadtxt(file, usecols=2))
    Cv.append(np.loadtxt(file, usecols=3))
    X.append(np.loadtxt(file, usecols=4))


maxX = np.argmax(X, axis=1)
m, b, r_value, p_value, std_err = stats.linregress(1./L,t[maxX])

fig = plt.figure()
plt.plot(1. / L, t[maxX], "o")
plt.plot(1. / L, m*(1. / L) + b)
plt.xlabel("$1/L$")
plt.ylabel("$T_c$ $[\,J/k_B\,]$")
plt.legend(["Finite Lattice", "%.4f x + %.4f $\\pm$ %s" % (m, b,std_err)])
plt.grid()
fig.savefig("plots/critTemp.pdf")

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, E[i], linewidth=0.8)
plt.xlabel("T $[\,J/k_B\,]$")
plt.ylabel("$\\langle E \\rangle /L^2$ $[\,J\,]$")
plt.legend(["40x40 Lattice", "60x60 Lattice", "100x100 Lattice",
            "160x160 Lattice"])
plt.grid()
fig.savefig("plots/evolution_energy.pdf")

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, M[i], linewidth=0.8)
plt.xlabel("T $[\,J/k_B\,]$")
plt.ylabel("$\\langle |M| \\rangle /L^2$")
plt.legend(["40x40 Lattice", "60x60 Lattice", "100x100 Lattice",
            "160x160 Lattice"])
plt.grid()
fig.savefig("plots/evolution_magnetization.pdf")

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, Cv[i], linewidth=0.8)
plt.xlabel("T $[\,J/k_B\,]$")
plt.ylabel("$ C_V /L^2$ $[\,k_B\,]$")
plt.legend(["40x40 Lattice", "60x60 Lattice", "100x100 Lattice",
            "160x160 Lattice"])
plt.grid()
fig.savefig("plots/evolution_cv.pdf")

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, X[i], linewidth=0.8)
plt.xlabel("T $[\,J/k_B\,]$")
plt.ylabel("$\\chi /L^2$")
plt.legend(["40x40 Lattice", "60x60 Lattice", "100x100 Lattice",
            "160x160 Lattice"])
plt.grid()
fig.savefig("plots/evolution_susceptibility.pdf")
# -----------------------------------------------------------------------------
