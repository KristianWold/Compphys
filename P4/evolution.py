import os
import numpy as np
import matplotlib.pyplot as plt
import math as m
import sys

cycles = int(sys.argv[1])
cutoff = int(sys.argv[2])
L = int(sys.argv[3])
T_start = float(sys.argv[4])
T_end = float(sys.argv[5])
T_step = float(sys.argv[6])

N = int(m.ceil((T_end - T_start) / T_step))

file = open("results/evolution_L=%s.txt" % L, "w")

for i in range(N + 1):
    os.system("mpirun -np 8 ./simulation.x %s %s %s %s" %
              (cycles, cutoff, L, (T_start + i * T_step)))
    expect = np.loadtxt("results/expection.txt", usecols=0)
    for val in expect:
        file.write("%s " % val)
    file.write("\n")
    print("go!")

file.close()
