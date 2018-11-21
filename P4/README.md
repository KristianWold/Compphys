# Project 4

This repository contains programs made for project 4 produced in a collaboration between [Lasse](https://github.com/lasselb87), [Nicolai](https://github.com/nicolossus), and [Kristian](https://github.com/KristianWold).


### Contents

The [Classes](https://github.com/KristianWold/Compphys/tree/master/P4/classes) folder contains object-oriented code
for calculating Ising Model estimations using the Metropolis Algorithm, with some accompaying unit tests.

The [Results](https://github.com/KristianWold/Compphys/tree/master/P4/results) and
[Plots](https://github.com/KristianWold/Compphys/tree/master/P4/plots)
folder contains all the results produced in this project.

**Usage**

To compile the c++ implementation of the Metropolis Algorithm, run `make simulation.x`.

To perform benchmarking, run `python evolution N 1000 2 1 3 0.01`, where N is the desired
number of cycles. Run then `python analyticalVsNumerical N` to produce the plot.

To produce the plots showing how the system moves toward equilibrium, run `python equilibrium.py`.

To produce the probability distribution, run `python distribution.py`.

To produce plots of the number of accepted states as a function of N and T, run `python acceptedStates.py`.

To produce the plots showing the phase transition and finite size scaling, run `python phaseTrans.py`.

To test some of the code, run `make test.x`.
