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
    double dt;

    void totalAcceleration(mat &totalAcc, vec acc(vec, vec));
    void coordinatesToFile(ofstream &myfile);
    void euler(vec acc(vec,vec));
    void verlet(vec acc(vec,vec));

    mat totalAcc;
    mat prevAcc;
    mat prevVel;

public:
    mat kineticEnergy;
    mat potentialEnergy;
    vec energyAllPlanets;
    vec angularMomentum;

    Solver(vector<Planet> p, double scale);

    void solve(int method, vec acc(vec, vec), double T, int N, int sampleN);

    void sampleEnergyAndAngular(mat &kinetic, mat &potential, vec &angular,
                                vec &energyAllPlanets, int i);

    double kineticFluctuation(int i);
    double potentialFluctuation(int i);
    double totalEnergyFluctuation();
    double angularFluctuation();
};

#endif
