import os

T_start = 2.1
T_end = 2.5
T_step = 0.01
N = int((T_end - T_start) / T_step)
file = open("results/phase.txt", "w")
file.close()

for i in range(N):
    os.system("mpirun -np 8 ./main.x 5000 100 %s" % (T_start + i * T_step))

os.system("python plot.py")
