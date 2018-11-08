import os

T_start = 2.1
T_end = 2.5
T_step = 0.05
N = int((T_end - T_start) / T_step)
file = open("results/phase.txt", "w")
file.close()

for i in range(N):
    os.system("mpirun -np 8 ./main.x 10000 100 %s" % (T_start + i * T_step))
    print("go!")

os.system("python plot.py")
