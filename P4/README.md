# Project 4

This repository contains programs made for project 4 produced in a collaboration between [Lasse](https://github.com/lasselb87), [Nicolai](https://github.com/nicolossus), and [Kristian](https://github.com/KristianWold).

In this project we simulate phase transitions with the Ising model in two dimensions. We use the Metropolis algorithm , a Markov chain Monte Carlo (MCMC) method, to calculate Ising model estimations.


### Contents

The [Classes](https://github.com/KristianWold/Compphys/tree/master/P4/classes) folder contains object-oriented code
for calculating Ising Model estimations using the Metropolis Algorithm. The folder contains unit tests verifying the object-oriented code as well. Files in the folder:

- [ising.hpp](https://github.com/KristianWold/Compphys/tree/master/P4/classes/ising.hpp): class framework for the Ising model.

- [ising.cpp](https://github.com/KristianWold/Compphys/tree/master/P4/classes/ising.cpp): member functions of the Ising class.

- [catch.hpp](https://github.com/KristianWold/Compphys/tree/master/P4/classes/catch.hpp): test framework for unit-tests, retrieved from the GitHub repository [catchorg/Catch2](https://github.com/catchorg/Catch2).

- [test.cpp](https://github.com/KristianWold/Compphys/tree/master/P4/classes/test.cpp): unit tests for verifying the Ising class.

- [metropolis.hpp](https://github.com/KristianWold/Compphys/tree/master/P4/classes/metropolis.hpp): class framework for the Metropolis algorithm.

- [metropolis.cpp](https://github.com/KristianWold/Compphys/tree/master/P4/classes/metropolis.cpp): member functions of the Metropolis class.

The [Results](https://github.com/KristianWold/Compphys/tree/master/P4/results) folder contains raw data results.

The [Plots](https://github.com/KristianWold/Compphys/tree/master/P4/plots) folder contains all figures produced in this project.

[simulation.cpp](https://github.com/KristianWold/Compphys/tree/master/P4/simulation.cpp): simulate Ising model estimations using the Metropolis algorithm. The code is parallelized with MPI.

[evolution.py](https://github.com/KristianWold/Compphys/tree/master/P4/evolution.py): a general recipe for running the simulation and writing the results to file for a range of temperatures.

[analyticalVsNumerical.py](https://github.com/KristianWold/Compphys/tree/master/P4/analyticalVsNumerical.py): benchmark simulated estimations with analytical results of a 2x2 lattice Ising model.

[equilibrium.py](https://github.com/KristianWold/Compphys/tree/master/P4/equilibrium.py): plotting procedures for plotting the equilibrium situation.

[distribution.py](https://github.com/KristianWold/Compphys/tree/master/P4/distribution.py): count relative frequency of unique energy states and plots the distribution.

[acceptedStates.py](https://github.com/KristianWold/Compphys/tree/master/P4/acceptedStates.py): plotting procedures for plotting the number of accepted states.

[phaseTrans.py](https://github.com/KristianWold/Compphys/tree/master/P4/phaseTrans.py): plotting procedures for plotting the thermodynamic quantities around a phase transition.


**Usage**

To compile the C++ implementation of the Metropolis Algorithm, run `make simulation.x`.

To perform benchmarking, run `python evolution.py N 1000 2 1 3 0.01`, where N is the desired
number of MC cycles. Run then `python analyticalVsNumerical.py N` to produce the plot.

To produce the plots showing how the system moves toward equilibrium, run `python equilibrium.py`.

To produce the probability distribution, run `python distribution.py`.

To produce plots of the number of accepted states as a function of N and T, run `python acceptedStates.py`.

To produce the plots showing the phase transition and finite size scaling, run `python phaseTrans.py`.

To test some of the features of the Ising class, run `make test.x`.
