#ifndef MONTECARLO_H
#define MONTECARLO_H

#include<stdlib.h>
#include<armadillo>
#include<fstream>
#include<cmath>
#include<map>
#include<random>
using namespace std;
using namespace arma;

class MonteCarlo
{
private:
    Spins spins;
    double acceptAmp;

public:
    int* energyAndMag;
    int E;
    int M;
    int accepted;

    uniform_real_distribution<float> rand_float;

    MonteCarlo(Spins spins);
    void solve(int cycles, mt19937 &engine);
};

#endif
