import numpy as np
import matplotlib.pyplot as plt
import struct

f = open("error.bin","rb")

#all generated vectors in one array
x = np.array([struct.unpack("d", f.read(8))[0] for i in range(7)])
y = np.array([struct.unpack("d", f.read(8))[0] for i in range(7)])

plt.plot(np.log10(x), np.log10(y))
#plt.plot(np.log10(x), np.log10(y))
plt.plot(np.log10(x), 2*np.log10(x))
plt.xlabel("log(h)")
plt.ylabel("log(relative error)")
plt.show()
