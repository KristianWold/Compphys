#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Set fontsizes in figures
params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
plt.rcParams.update(params)
plt.gcf().set_tight_layout(True)

# Benchmarking for non-interacting case. Different alpha
# ---------------------------------------------------------------------------
filename = "results/noninteracting.txt"
alpha = np.loadtxt(filename, usecols=0)
x = np.linspace(alpha[0], alpha[-1], 1000)
E = np.loadtxt(filename, usecols=1)
Var = np.loadtxt(filename, usecols=2)

fig = plt.figure(figsize=(10, 8))
plt.subplot(2, 1, 1)
plt.plot(alpha, E)
plt.plot(x, 3. / (2 * x) * (1 - x**2) + 3 * x, "--")
plt.xlabel("$\\alpha$")
plt.ylabel("$\\langle E \\rangle$")
plt.legend(["MC simulation", "Analytical"])
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(alpha, Var)
plt.xlabel("$\\alpha$")
plt.ylabel("$\\sigma^2$")
plt.legend(["MC simulation"])
plt.grid()

fig.savefig("./plots/noninteracting.pdf")
# ---------------------------------------------------------------------------

# Benchmarking for interacting case. Different alpha
# ---------------------------------------------------------------------------
filename = "results/interacting.txt"
alpha = np.loadtxt(filename, usecols=0)
E = np.loadtxt(filename, usecols=1)
Var = np.loadtxt(filename, usecols=2)

fig = plt.figure(figsize=(10, 8))
plt.subplot(2, 1, 1)
plt.plot(alpha, E)
plt.xlabel("$\\alpha$")
plt.ylabel("$\\langle E \\rangle$")
plt.legend(["MC simulation"])
plt.grid()

plt.subplot(2, 1, 2)
plt.plot(alpha, Var)
plt.xlabel("$\\alpha$")
plt.ylabel("$\\sigma^2$")
plt.legend(["MC simulation"])
plt.grid()

fig.savefig("./plots/interacting.pdf")
# ---------------------------------------------------------------------------

# Virial theorem for interacting and non-interacting case.
# ---------------------------------------------------------------------------
filename = "results/virial.txt"
omega = np.loadtxt(filename, usecols=0)
ratio1 = np.loadtxt(filename, usecols=1)
ratio2 = np.loadtxt(filename, usecols=2)

fig = plt.figure(figsize=(10, 8))

plt.plot(omega, ratio1)
plt.plot(omega, ratio2)
plt.xlabel("$\\omega")
plt.ylabel("$\\langle T \\rangle/\\langle V \\rangle$")
plt.legend(["Non-interacting", "Interacting"])
plt.grid()

fig.savefig("./plots/virial.pdf")
# ---------------------------------------------------------------------------

# Electron cloud
# ---------------------------------------------------------------------------
filename = "results/data.dat"
x = np.loadtxt(filename, usecols=0)
y = np.loadtxt(filename, usecols=1)
z = np.loadtxt(filename, usecols=2)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

plt.plot(x, y, z, "o", alpha=0.3, markersize=0.3)
plt.xlabel("$\\omega")
plt.ylabel("$\\langle T \\rangle/\\langle V \\rangle$")
plt.legend(["Non-interacting", "Interacting"])
plt.grid()

fig.savefig("./plots/electronCloud.pdf")
# ---------------------------------------------------------------------------
