import numpy as np
import matplotlib.pyplot as plt

file = "data.txt"
y = np.loadtxt(file, usecols=0)
plt.plot(y)


#file = "magnetization.txt"
#t = np.loadtxt(file, usecols=0)
#M = np.loadtxt(file, usecols=1)
#plt.plot(t, M)

plt.show()
