import os
import numpy as np
import matplotlib.pyplot as plt
import math as m
import sys

L = int(sys.argv[1])

T_start = 2.2
T_end = 2.4
T_step = 0.01
N = int(m.ceil((T_end - T_start) / T_step))

file = open("results/evolution_L=%s.txt"%L,"w")

for i in range(N+1):
    os.system("mpirun -np 8 ./simulation.x 1000000 1000 %s %s" %(L, (T_start + i * T_step)))
    expect = np.loadtxt("results/expection.txt",usecols=0)
    for val in expect:
        file.write("%s "%val)
    file.write("\n")
    print("go!")

file.close()
