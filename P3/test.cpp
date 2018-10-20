#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include <iostream>
#include <armadillo>
#include <fstream>
#include <vector>
#include "Classes/solver.h"
#include "Classes/planet.h"
using namespace std;
using namespace arma;

double scale = 4*M_PI*M_PI/(365.25*365.25);

vec newton(vec pos, vec vel)
{
    double rCube = pow(norm(pos), 3);
    return -scale/rCube*pos;
}

TEST_CASE("Circular orbit")
//Check that the kinetic and potential energy, and angular momentum as well,
//are conserved independently for circular orbit
{
    double eps = 1e-5;

    Planet Earth(vec({1,0,0}), vec({0,2*M_PI/365.25,0}),1./333000);
    vector<Planet> solarsystem = vector<Planet>{Earth};
    Solver solver(solarsystem, scale);

    solver.solve(2, newton, 365.25, 1000000, 1000);

    REQUIRE(solver.kineticFluctuation(0) < eps);
    REQUIRE(solver.potentialFluctuation(0) < eps);
    REQUIRE(solver.angularFluctuation() < eps);
}

TEST_CASE("Elliptic orbit")
//Check that combined kinetic and potential energy is conserved
{
    double eps = 1e-5;

    Planet Earth(vec({1,0,0}), vec({0,M_PI/365.25,0}),1./333000);
    vector<Planet> solarsystem = vector<Planet>{Earth};
    Solver solver(solarsystem, scale);

    solver.solve(2 ,newton, 365.25, 1000000, 1000);

    REQUIRE(solver.totalEnergyFluctuation() < eps);
    REQUIRE(solver.angularFluctuation() < eps);
}

TEST_CASE("Three-body")
//Check that combined kinetic and potential energy between both planets
//are conserved
{
    double eps = 1e-5;

    Planet Earth(vec({1,0,0}), vec({0,2*M_PI/365.25,0}),1./333000);
    Planet Jupiter(vec({3,0,0}), vec({0,M_PI/365.25,0}),0.09);
    vector<Planet> solarsystem = vector<Planet>{Earth,Jupiter};
    Solver solver(solarsystem, scale);

    solver.solve(2, newton, 365.25, 1000000, 1000);
    REQUIRE(solver.totalEnergyFluctuation() < eps);
    REQUIRE(solver.angularFluctuation() < eps);
}
