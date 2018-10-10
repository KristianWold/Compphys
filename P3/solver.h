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

class Verlet
{
private:
    //keeps track of if the methods solve() and solveEnergy() have been called
    bool solved = false;
    bool solvedEnergy = false;

    vec t;
    cube pos;
    cube vel;
    vector<Planet> planets;
    int numPlanets;
    int N;
    double scale;

    void totalAcceleration(mat &totalAcc, vec acc(vec, vec), int i);

public:
    mat kineticEnergy;
    mat potentialEnergy;
    mat angularMomentum;

    Verlet(vector<Planet> p, int n, double scale);

    void solve(vec acc(vec, vec), double T, double dt);

    void solveEnergy();

    void coordinatesToFile(string filename);

    void energyToFile(string filename);
};

#endif
