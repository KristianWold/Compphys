import numpy as np
import matplotlib.pyplot as plt
import sys

file = "data.txt"

x = np.loadtxt(file, usecols=0)
y = np.loadtxt(file, usecols=1)
z = np.loadtxt(file, usecols=2)

r = np.sqrt(x*x + y*y + z*z)

atPerihelion = False
rPerihelion = 0.3075
eps = 1e-4

for i in range(len(x)):
    if abs((r[i]-rPerihelion)/rPerihelion) < eps and not atPerihelion:
        atPerihelion = True
        angle.append(y[i]/x[i])

    if abs((r[i]-rPerihelion)/rPerihelion) > eps and atPerihelion:
        atPerihelion = False

print np.arctan(angle)
