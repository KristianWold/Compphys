import numpy as np
import matplotlib.pyplot as plt

meta = "results/meta.txt"
metainfo = np.loadtxt(meta, usecols=0)

cycles = int(metainfo[0])
cores = int(metainfo[1])
L = int(metainfo[2])
T = float(metainfo[3])

file = "results/data.dat"
array = np.fromfile(file, dtype="int32", count=-1)

E = np.zeros(cores * cycles)
M = np.zeros(cores * cycles)

for i in range(cores):
    start = 2 * cycles * i
    end = 2 * cycles * i + cycles
    E[(i * cycles): ((i + 1) * cycles)] = array[start:end]
    M[(i * cycles): ((i + 1) * cycles)] = array[(start + cycles):(end + cycles)]

E_mean = sum(E) / (cores * cycles)
M_mean = sum(M) / (cores * cycles)

Cv = 1 / T**2 * (sum(E**2) / (cores * cycles) - E_mean**2)
X = 1 / T * (sum(M**2) / (cores * cycles) - M_mean**2)

print("Energy = %s" % E_mean)
print("Magnetization = %s" % M_mean)
print("Heat Capacity %s" % Cv)
print("Suceptibility %s" % X)

"""
P = {}

for i in range(cycles * cores):
    energy = E[i]
    if energy in P:
        P[energy] += 1
    else:
        P[energy] = 1

plt.plot(P.keys(), P.values(), "o")
plt.show()
"""
