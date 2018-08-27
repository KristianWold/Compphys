import numpy as np
import matplotlib.pyplot as plt
import struct

f = open("error.bin","rb")

#all generated vectors in one array
x = [struct.unpack("d", f.read(8)) for i in range(30)]
y = [struct.unpack("d", f.read(8)) for i in range(30)]

#analytical solution
plt.plot(x, y)
plt.show()
