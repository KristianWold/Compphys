import numpy as np
import matplotlib.pyplot as plt

file = "data.dat"
y = np.fromfile(file, dtype="int32", count=-1)
plt.plot(y)
plt.show()

cutoff = 0

print(sum(y[cutoff:-1]) / float(len(y) - cutoff))
