import numpy as np
import matplotlib.pyplot as plt
"""
meta = "results/meta.txt"
metainfo = np.loadtxt(meta, usecols=0)

cycles = int(metainfo[0])
cores = int(metainfo[1])
L = int(metainfo[2])
T = float(metainfo[3])

file = "results/data.dat"
array = np.fromfile(file, dtype="int32", count=-1)

for i in range(cores):
    start = int(2 * (i * cycles))
    end = int(2 * (i * cycles) + cycles)
    plt.figure()
    plt.plot(array[start: end])
"""
file = "results/phase.txt"
T = np.loadtxt(file, usecols=0)
M = np.loadtxt(file, usecols=1)
plt.plot(T, M)
plt.xlabel("T")
plt.ylabel("M")
plt.show()
