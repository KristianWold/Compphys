import numpy as np
import matplotlib.pyplot as plt
import struct

meta = open("meta.bin","rb")
f = open("myfile.bin","rb")

#number of vectors
N = struct.unpack("i", meta.read(4))[0]

#dimmension of vectors
n = [struct.unpack("i", meta.read(4))[0] for i in range(N)]

#all generated vectors in one array
y = [struct.unpack("d", f.read(8)) for i in range(sum(n))]

#analytical solution
x = np.linspace(0,1,1000)
def solution(x):
    return 1 - (1-np.exp(-10))*x -np.exp(-10*x)

x0 = 0
for i in range(N):
    #plots a slice of y corresponding to the correct vector
    fig = plt.figure("n%s"%n[i])
    ax = fig.gca()
    plt.grid()
    plt.plot(np.linspace(0,1,n[i]), y[x0:(n[i]+x0)])
    plt.plot(x, solution(x))
    plt.legend(["n=%s"%n[i], "analytical"])


    fig.savefig("results/n%s"%n[i])
    #plt.show()
    #Advances the slice to the next vector
    x0 += n[i]
