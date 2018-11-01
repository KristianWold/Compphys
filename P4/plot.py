import numpy as np
import matplotlib.pyplot as plt

file = "data.dat"

N = 100000

EM = np.fromfile(file, dtype="int32", count=-1)

plt.plot(EM[N:-1])
plt.show()

cutoff = 3000

print(sum(E[cutoff:N]) / float(N - cutoff))
