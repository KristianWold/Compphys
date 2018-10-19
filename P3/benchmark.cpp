#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include "solver.h"
#include "planet.h"
using namespace std;
using namespace arma;

double scale = 4*M_PI*M_PI/(365.25*365.25);

vec newton(vec pos, vec vel)
{
    double rCube = pow(norm(pos), 3);
    return -scale/rCube*pos;
}

int main()
{
    Planet Earth(vec({1,0,0}), vec({0,2*M_PI/365.25,0}),1./333000);
    vector<Planet> solarsystem = vector<Planet>{Earth,Jupiter};

    double T = 1000; 



    return 0;
}
