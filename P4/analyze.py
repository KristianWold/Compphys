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

cutoff = [10, 10, 10, 10, 10, 10, 10, 10]

E = 0
M = 0
E2 = 0
M2 = 0
normalize = 0.

for i in range(cores):
    start = int(2 * (i * cycles) + cutoff[i])
    end = int(2 * (i * cycles) + cycles)
    E += sum(array[start:end])
    E2 += sum(array[start:end]**2)
    M += sum(abs(array[(start + cycles):(end + cycles)]))
    M2 += sum(array[(start + cycles):(end + cycles)]**2)
    normalize += (cycles - cutoff[i])

E /= normalize
M /= normalize
E2 /= normalize
M2 /= normalize


Cv = np.sqrt(E2 - E**2)
X = np.sqrt(M2 - M**2)

print("Mean E=%s" % E)
print("Mean M=%s" % M)
