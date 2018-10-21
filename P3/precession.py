import numpy as np
import matplotlib.pyplot as plt
import sys

file = "data.txt"

x = np.loadtxt(file, usecols=0)
y = np.loadtxt(file, usecols=1)
z = np.loadtxt(file, usecols=2)

r = np.sqrt(x*x + y*y + z*z)

angle = []

for i in range(1,len(x)-1):
    if ((r[i-1] > r[i]) and (r[i+1] > r[i])):
        angle.append(y[i]/x[i])


#plt.plot(x,y)

plt.plot(np.arctan(angle))
print(len(angle))
plt.show()
