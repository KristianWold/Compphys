import os
import numpy as np
import matplotlib.pyplot as plt
import math as m
import sys

# phase transtion
# -----------------------------------------------------------------------------
L = [40, 60, 80, 100, 140, 160]
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

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, E[i])
plt.xlabel("T")
plt.ylabel("<E>")
plt.legend(["L=40", "L=60", "L=80", "L=100", "L=140"])
fig.savefig("plots/evolution_energy.pdf")

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, M[i])
plt.xlabel("temperature")
plt.ylabel("<|M|>")
plt.legend(["L=40", "L=60", "L=80", "L=100", "L=140"])
fig.savefig("plots/evolution_magnetization.pdf")

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, Cv[i])
plt.xlabel("Temperature")
plt.ylabel("Cv")
plt.legend(["L=40", "L=60", "L=80", "L=100", "L=140"])
fig.savefig("plots/evolution_cv.pdf")

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, X[i])
plt.xlabel("temperature")
plt.ylabel("Susceptibility")
plt.legend(["L=40", "L=60", "L=80", "L=100", "L=140"])
fig.savefig("plots/evolution_susceptibility.pdf")
# -----------------------------------------------------------------------------
