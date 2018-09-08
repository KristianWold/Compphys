import numpy as np
import matplotlib.pyplot as plt
import struct

f = open("error.bin","rb")

#different stepsize
h = np.array([struct.unpack("d", f.read(8))[0] for i in range(7)])

#maximum error for each stepsize
err = np.array([struct.unpack("d", f.read(8))[0] for i in range(7)])

fig = plt.figure()
ax = fig.gca()
plt.plot(np.log10(h), np.log10(err))
#graph with slope of 2
plt.plot(np.log10(h), 2*np.log10(h))
plt.xlabel("log(h)")
plt.ylabel("log(relative error)")
plt.legend(["Relative error of numerical method","y = 2x"])
fig.savefig("results/errorestimate.png")
plt.show()
