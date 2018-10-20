import numpy as np
import matplotlib.pyplot as plt
import sys

# Set fontsizes in figures
params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large'}
plt.rcParams.update(params)

if sys.argv[1] == "pos":
    file = "data.txt"

    x1 = np.loadtxt(file, usecols=0)
    y1 = np.loadtxt(file, usecols=1)

    #x2 = np.loadtxt(file, usecols=3)
    #y2 = np.loadtxt(file, usecols=4)

    fig = plt.figure()
    plt.gca().set_aspect("equal")
    plt.plot(x1, y1)
    #plt.plot(x2, y2)
    plt.show()

if sys.argv[1] == "energy":
    file = "energy.txt"

    x = np.loadtxt(file, usecols=0)
    y = np.loadtxt(file, usecols=1)

    fig = plt.figure()
    #plt.gca().set_aspect("equal")
    plt.plot(np.linspace(0,1,len(x)),x)
    plt.plot(np.linspace(0,1,len(y)),y)
    plt.grid(True)
    plt.show()

if sys.argv[1] == "fluctuation":
    file = "fluctuation.txt"

    n = np.loadtxt(file, usecols=0)
    y = np.loadtxt(file, usecols=1)

    fig = plt.figure()
    plt.plot(np.log10(n),np.log10(y))
    plt.grid(True)
    plt.show()


"""
x1 = np.loadtxt(file, usecols=0)
y1 = np.loadtxt(file, usecols=1)

x2 = np.loadtxt(file, usecols=6)
y2 = np.loadtxt(file, usecols=7)

#x3 = np.loadtxt(file, usecols=12)
#y3 = np.loadtxt(file, usecols=13)



plt.plot(x1, y1)
plt.plot(x2, y2)
#plt.plot(x3, y3)
plt.show()
"""
