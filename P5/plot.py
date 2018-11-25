import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")

x1 = np.loadtxt("data.txt", usecols=(0, 1, 2))
x2 = np.loadtxt("data.txt", usecols=(3, 4, 5))

a = np.loadtxt("moredata.txt", usecols=0)
E1 = np.loadtxt("moredata.txt", usecols=1)
E2 = np.loadtxt("moredata.txt", usecols=2)
E3 = np.loadtxt("moredata.txt", usecols=3)

# plt.plot(x1[:, 0], x1[:, 1], x1[:, 2], 'o', c="b", alpha=0.01, markersize=1)
# plt.plot(x2[:, 0], x2[:, 1], x2[:, 2], 'o', c="r", alpha=0.01, markersize=1)
plt.plot(a, E1)
plt.plot(a, E2)
plt.plot(a, E3)
plt.plot(a, 3. / (2 * a) * (1 - a**2) + 3 * a)
plt.xlabel("$\\alpha$")
plt.ylabel("$\\langle E \\rangle$")
plt.legend(["Non-interacting", "Interacting", "Pade-Jastrow", "Analytical"])
plt.show()
