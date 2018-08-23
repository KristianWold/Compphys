#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def velfield(n=20):
    x = np.linspace(0, 2 * np.pi, n)
    y = np.linspace(0, 2 * np.pi, n)
    [X, Y] = np.meshgrid(x, y)
    u = np.cos(X) * np.sin(Y)
    v = -np.sin(X) * np.cos(Y)
    return x, y, u, v


n = 1000
x, y, u, v = velfield(n)
plt.axis('off')
plt.plot(u, v)
plt.show()
