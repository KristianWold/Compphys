import os
import numpy as np
import matplotlib.pyplot as plt
import math as m
import sys

accepted = []
N = [1e3, 1e4, 1e5, 1e6, 1e7]

for n in N:
    os.system("mpirun -np 8 ./simulation.x %s 0 20 2" % n)
    accepted.append(np.loadtxt("results/meta.txt", usecols=0)[5])

fig = plt.figure()
plt.plot(N, accepted)
plt.xlabel("Cycles")
plt.ylabel("Accepted States")
plt.legend("Accepted states, L=20, T = 2")
fig.savefig("plots/acceptedCycles.pdf")

accepted = []
T = np.linspace(1, 3, 20)
for t in T:
    os.system("mpirun -np 8 ./simulation.x 10000 0 20 %s" % t)
    accepted.append(np.loadtxt("results/meta.txt", usecols=0)[5])

fig = plt.figure()
plt.plot(T, accepted)
plt.xlabel("T")
plt.ylabel("Accepted States")
plt.legend("Accepted states, L=20, cycles = 10000")
fig.savefig("plots/acceptedTemp.pdf")
