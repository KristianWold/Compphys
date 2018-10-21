import numpy as np
import matplotlib.pyplot as plt
import sys

file1 = "angle1.txt"

angle = np.loadtxt(file1, usecols=0)

fig = plt.figure()
plt.plot(np.arctan(angle))
fig.savefig("results/precessionNewton.pdf")

file2 = "angle2.txt"
angle = np.loadtxt(file2, usecols=0)

fig = plt.figure()
plt.plot(np.arctan(angle))
fig.savefig("results/precessionEinstein.pdf")
