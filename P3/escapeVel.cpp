#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include "Classes/solver.h"
#include "Classes/planet.h"
using namespace std;
using namespace arma;

double scale = 4*M_PI*M_PI;

vec newton(vec pos, vec vel)
{
    double rCube = pow(norm(pos), 3);
    return -scale/rCube*pos;
}

int main(int argc, char const *argv[])
{
    Planet Earth(vec({1,0,0}), vec({0,sqrt(2)*2*M_PI,0}),1./333000);
    vector<Planet> solarsystem = vector<Planet>{Earth};
    Solver solver(solarsystem, scale);

    double T = atof(argv[1]);
    int N = atoi(argv[2]);
    int sampleN = atoi(argv[3]);

    solver.solve(2, newton, T, N, sampleN);
    system("python3 plot.py pos");

    return 0;
}