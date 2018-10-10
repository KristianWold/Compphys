#include "planet.h"
using namespace std;
using namespace arma;

Planet::Planet()
{
    pos = zeros(3);
    vel = zeros(3);
    M = 0;
}
Planet::Planet(vec position, vec velocity, double mass)
{
    pos = position;
    vel = velocity;
    M = mass;
}
