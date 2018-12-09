#ifndef VMC_HPP
#define VMC_HPP

#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <armadillo>
#include <fstream>

using namespace std;
using namespace chrono;
using namespace arma;

//function pointer type to pass special functons to the class VMC
typedef double (*myfunc1)(double*, double, mat&, vec&, int);
typedef double (*myfunc2)(double*, double, mat&);

struct Result
{
    double E;
    double Var;
    double kinetic;
    double potential;
    double R12;
    int accepted;
};

class VMC
{
public:
    myfunc1 acceptAmp;
    myfunc2 localKinetic, localPotential;
    int numDim, numParticles, accepted;
    int numCycles, preCycles;
    double kineticE, potentialE, k_E, p_E, E, E2, Var, R12, step;
    double* params;
    double omega;

    mat pos;
    vec delta;

    mt19937_64 engine;
    uniform_real_distribution<double> myRandu = uniform_real_distribution<double>(0,1);

    VMC(int numDim, int numParticles, myfunc1 acceptAmp, myfunc2 localKinetic, myfunc2 localPotential);

    void mcCycle();

    Result solve(int numCycles, int preCycles, double* params, double omega, bool writeToFile);

    void optimize(double *params, int numParams, double range, int step,
                  int maxIter, int numCycles, int preCycles);
};

#endif
