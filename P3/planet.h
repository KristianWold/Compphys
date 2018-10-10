#ifndef PLANET_H
#define PLANET_H

#include <armadillo>
using namespace std;
using namespace arma;

class Planet
{
public:
    vec pos;
    vec vel;
    double M;
    Planet();
    Planet(vec position, vec velocity, double mass);
};

#endif
