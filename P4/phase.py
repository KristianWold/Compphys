import os
import numpy as np
import matplotlib.pyplot as plt

T_start = 2.24
T_end = 2.3
T_step = 0.01
N = int((T_end - T_start) / T_step)

T = []
X = []

for i in range(N+1):
    os.system("mpirun -np 8 ./simulation.x 30000 140 %s" % (T_start + i * T_step))
    expect = np.loadtxt("results/expection.txt", usecols=0)
    T.append(expect[0])
    X.append(expect[3])
    print("go!")

np.savetxt("results/phase.txt", np.transpose((T, X)))

plt.plot(T,X)
plt.show()
