import numpy as np
import matplotlib.pyplot as plt
import sys

file1 = "angle1.txt"

radians = np.arctan(np.loadtxt(file1, usecols=0))
arcseconds = radians * 180 * 3600 / np.pi
t = np.linspace(0, 100, len(arcseconds))


fig = plt.figure()
fit = np.poly1d(np.polyfit(t, arcseconds, 1))
plt.plot(t, arcseconds)
plt.plot(t, fit(t))
plt.xlabel("Year")
plt.xlabel("arcseconds")
plt.legend(["precession", fit])
fig.savefig("results/precessionNewton.pdf")


file2 = "angle2.txt"
radians = np.arctan(np.loadtxt(file2, usecols=0))
arcseconds = radians * 180 * 3600 / np.pi
t = np.linspace(0, 100, len(arcseconds))

fig = plt.figure()
fit = np.poly1d(np.polyfit(t, arcseconds, 1))
plt.plot(t, arcseconds)
plt.plot(t, fit(t))
plt.xlabel("Year")
plt.ylabel("arcseconds")
plt.legend(["precession", fit])
fig.savefig("results/precessionEinstein.pdf")

precession = fit(100) - fit(0)
print(precession)
