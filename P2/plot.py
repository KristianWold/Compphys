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

if sys.argv[1] == "benchmark":
    file = "benchmark.txt"

    n = np.log10(np.loadtxt(file, usecols=0))
    t1 = np.log10(np.loadtxt(file, usecols=1))
    t2 = np.log10(np.loadtxt(file, usecols=2))
    N = np.log10(np.loadtxt(file, usecols=3))

    fit1 = np.poly1d(np.polyfit(n, t1, 1))
    fit2 = np.poly1d(np.polyfit(n, t2, 1))
    fit3 = np.poly1d(np.polyfit(n, N, 1))

    fig1 = plt.figure(1)
    title = "Jacobi's Method vs Armadillo's eig_sym()"

    plt.plot(n, t1)
    plt.plot(n, fit1(n), ":")
    plt.plot(n, t2)
    plt.plot(n, fit2(n), ":")

    plt.gca().set_xlabel('$\\log_{10}(n)$')
    plt.gca().set_ylabel('$\\log_{10}(t)$')
    plt.gca().set_title(title)
    plt.legend(["Jacobi", fit1, "Arma", fit2])
    fig1.savefig("./results/benchmark.png")

    fig2 = plt.figure(2)
    title = "Number of iterations for Jacobi's method"
    plt.plot(n, N)
    plt.plot(n, fit3(n), ":")

    plt.gca().set_xlabel('$\\log_{10}(n)$')
    plt.gca().set_ylabel('$\\log_{10}(N)$')
    plt.gca().set_title(title)
    plt.legend(["Number of iterations", fit3])
    fig2.savefig("./results/benchmark_iter.png")

if sys.argv[1] == "eigenvectors":
    file = "eigenvec.txt"

    w = sys.argv[2]
    x = np.loadtxt(file, usecols=0)
    y1 = np.loadtxt(file, usecols=1)
    y2 = np.loadtxt(file, usecols=2)

    fig = plt.figure(1)
    title = "Harmonic Oscillator Well with Two Electrons"

    plt.plot(x, y1**2)
    plt.plot(x, y2**2)

    plt.gca().set_xlabel('$\\rho$')
    plt.gca().set_ylabel('Probability')
    plt.gca().set_title(title)

    plt.gcf().set_tight_layout(True)
    x_pos = 0.8
    y_pos = 0.95
    plt.text(x_pos, y_pos, "$\\omega=$" + w,
         horizontalalignment='center',
         verticalalignment='center',
         transform=plt.gca().transAxes)


    plt.legend(["Non-interacting", "Interacting"])
    fig.savefig("./results/eingenvectors_w=%s.png" % w)
    plt.show()
