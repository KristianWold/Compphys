import numpy as np
import matplotlib.pyplot as plt
import sys

T = float(sys.argv[1])

E = -(8*np.sinh(8/T))/(np.cosh(8/T) + 3)
M = (2*np.exp(8/T) + 4)/(np.cosh(8/T) + 3)
Cv = 1/T**2*64*(3*np.cosh(8/T) + 1)/(np.cosh(8/T) + 3)**2
X = 1/T*(8*np.exp(8/T) + 8)/(np.cosh(8/T) + 3)

print("For T = %s"%T)
print("Energy = %s"%E)
print("Magnetization = %s"%M)
print("Cv = %s"%Cv)
print("X = %s"%X)
