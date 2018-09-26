import numpy as np
import matplotlib.pyplot as plt
import sys

if sys.argv[1] == "benchmark":
    file = "benchmark.txt"

    n = np.log10(np.loadtxt(file, usecols=0))
    t1 = np.log10(np.loadtxt(file, usecols=1))
    t2 = np.log10(np.loadtxt(file, usecols=2))

    fit1 = np.poly1d(np.polyfit(n, t1, 1))
    fit2 = np.poly1d(np.polyfit(n, t2, 1))

    fig = plt.figure("Jacobis method vs Armadillo eig_sym")
    plt.plot(n, t1)
    plt.plot(n, fit1(n), ":")
    plt.plot(n, t2)
    plt.plot(n, fit2(n), ":")
    plt.xlabel("log10(n)")
    plt.ylabel("log10(t)")
    plt.legend(["Jacobi", fit1, "Arma", fit2])
    fig.savefig("./results/benchmark.png")


if sys.argv[1] == "eigenvectors":
    file = "eigenvec.txt"

    w = sys.argv[2]
    x = np.loadtxt(file, usecols=0)
    y1 = np.loadtxt(file, usecols=1)
    y2 = np.loadtxt(file, usecols=2)

    fig = plt.figure("w = " + w)
    plt.plot(x, y1**2)
    plt.plot(x, y2**2)

    plt.xlabel("|x1 - x2|")
    plt.ylabel("probability")
    plt.legend(["non-interacting", "interacting"])
    fig.savefig("./results/eingenvectors_w=%s.png" % w)
    plt.show()
