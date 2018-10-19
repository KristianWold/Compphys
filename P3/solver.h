#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <armadillo>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include "planet.h"
using namespace std;
using namespace arma;

class Solver
{
private:
    //keeps track of if the methods solve() and solveEnergy() have been called
    bool solved = false;
    bool solvedEnergy = false;

    vec t;
    mat pos;
    mat vel;
    vector<Planet> planets;
    int numPlanets;
    int N;
    double scale;
    int countSample;

    void totalAcceleration(mat &totalAcc, vec acc(vec, vec));
    void coordinatesToFile(ofstream &myfile);

public:
    mat kineticEnergy;
    mat potentialEnergy;
    mat angularMomentum;

    Solver(vector<Planet> p, double scale);

    void solveVerlet(vec acc(vec, vec), double T, int N, int sampleN);

    void sampleEnergyAndAngular(mat &kinetic, mat &potential, mat &angular, int i);
};

#endif
