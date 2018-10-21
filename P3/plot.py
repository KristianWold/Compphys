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


    x = np.loadtxt(file, usecols=0)
    y = np.loadtxt(file, usecols=1)

    #x2 = np.loadtxt(file, usecols=3)
    #y2 = np.loadtxt(file, usecols=4)

    fig = plt.figure()
    plt.gca().set_aspect("equal")
    plt.plot(x, y)
    #plt.plot(x2, y2)
    plt.show()

if sys.argv[1] == "singlePlanet":
    file = "data.txt"

    name = sys.argv[2]
    T = sys.argv[3]
    N = sys.argv[4]

    x = np.loadtxt(file, usecols=0)
    y = np.loadtxt(file, usecols=1)

    fig = plt.figure()
    plt.gca().set_aspect("equal")
    plt.xlabel("AU")
    plt.ylabel("AU")
    plt.legend(["T=%s yr, N=%s, using %s"%(T,N,name)])
    plt.plot(x,y)
    fig.savefig("results/T=%s_N=%s_%s.pdf"%(T,N,name))

if sys.argv[1] == "energy":
    file = "energy.txt"

    T = float(sys.argv[2])
    N = int(sys.argv[3])

    x = np.loadtxt(file, usecols=0)
    y = np.loadtxt(file, usecols=1)

    fig = plt.figure()
    plt.plot(np.linspace(0,T,len(x)),x)
    plt.plot(np.linspace(0,T,len(x)),y)
    plt.grid(True)
    plt.xlabel("T")
    plt.ylabel("E")
    plt.legend(["Euler","Verlet"])
    fig.savefig("results/EulerVsVerlet_T=5_N=10000.pdf")

if sys.argv[1] == "fluctuation":
    file = "fluctuation.txt"

    name = sys.argv[2]

    n = np.loadtxt(file, usecols=0)
    y = np.loadtxt(file, usecols=1)

    fig = plt.figure()
    plt.plot(np.log10(n),np.log10(y))
    plt.xlabel("log10(N)")
    plt.ylabel("log10($\\epsilon$)")
    plt.grid(True)
    fig.savefig("results/fluctuation_%s.pdf"%name)


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
