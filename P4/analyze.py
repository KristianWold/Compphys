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

cutoff = 500
length = (cycles - cutoff) * cores

M = np.zeros((cores, cycles))
E = np.zeros((cores, cycles))


for i in range(cores):
    start = 2 * cycles * i
    end = 2 * cycles * i + cycles
    E[i] = array[start:end]
    M[i] = array[(start + cycles):(end + cycles)]

# plt.plot(E[0])
# plt.show()
E_mean = 0
M_mean = 0
E2 = 0
M2 = 0

for i in range(cores):
    truncE = E[i, cutoff:cycles]
    truncM = M[i, cutoff:cycles]

    E_mean += sum(truncE)
    M_mean += sum(truncM)

    E2 += sum(truncE**2)
    M2 += sum(truncM**2)


E_mean /= length * L * L
M_mean /= length * L * L
E2 /= length
M2 /= length

Cv = 1 / T**2 * (E2 - E_mean**2)
X = 1 / T * (M2 - M_mean**2)

print("T = %s" % T)
print("Energy = %s" % E_mean)
print("Magnetization = %s" % M_mean)
file = open("results/phase.txt", "a")
file.write("%s %s\n" % (T, M_mean))
file.close()

#print("Heat Capacity %s" % Cv)
#print("Suceptibility %s" % X)


P = {}

for i in range(cores):
    for j in range(cycles):
        energy = E[i, j]
        if energy in P:
            P[energy] += 1
        else:
            P[energy] = 1

for p in P:
    P[p] /= cycles * cores


#plt.plot(P.keys(), P.values(), "o")
# plt.show()
