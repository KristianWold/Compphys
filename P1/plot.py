import numpy as np
import matplotlib.pyplot as plt
import struct
import sys

metaName = sys.argv[1]
myfileName = sys.argv[2]

meta = open(metaName + ".bin","rb")
f = open(myfileName + ".bin","rb")

#number of vectors
N = struct.unpack("i", meta.read(4))[0]

#dimmension of vectors
n = [struct.unpack("i", meta.read(4))[0] for i in range(N)]

#all generated vectors in one array
y = [struct.unpack("d", f.read(8)) for i in range(sum(n))]

x0 = 0
for i in range(N):
    #plots a slice of y corresponding to the correct vector
    plt.plot(np.linspace(0,1,n[i]), y[x0:(n[i]+x0)])
    #print(y[x0:(n[i]+x0)])
    x0 += n[i]

#analytical solution
x = np.linspace(0,1,1000)
plt.plot(x, 1 - (1-np.exp(-10))*x -np.exp(-10*x))

plt.legend(["n=10","n=100","n=1000", "analytical"])
plt.show()
