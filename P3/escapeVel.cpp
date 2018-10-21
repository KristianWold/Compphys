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
    double T = 5;
    int N = 100000;
    int sampleN = 10;

    Planet Earth(vec({1,0,0}), vec({0,0,0}),1./333000);
    vector<Planet> solarsystem = vector<Planet>{Earth};

    for(int i = 0; i<6; i++)
    {
        solarsystem[0].vel(1) = sqrt(1+0.25*i)*2*M_PI;
        Solver solver(solarsystem, scale);
        solver.solve(2, newton, T, N, sampleN, "data" + to_string(i) + ".txt");
    }
    system("python3 plot.py escapeVel");
    return 0;
}
