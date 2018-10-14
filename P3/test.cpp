#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include <iostream>
#include <armadillo>
#include <fstream>
#include <vector>
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

TEST_CASE("Circular orbit")
//Check that the kinetic and potential energy, and angular momentum as well,
//are conserved independently for circular orbit
{
    double eps = 1e-5;

    Planet Earth(vec({1,0,0}), vec({0,2*M_PI/365.25,0}),1./333000);
    vector<Planet> solarsystem = vector<Planet>{Earth};
    Verlet solver(solarsystem, 1, scale);

    solver.solve(newton, 365.25, 1000000, 1000);

    vec kinetic = solver.kineticEnergy.col(0);
    double maxKinetic = kinetic(kinetic.index_max());
    double minKinetic = kinetic(kinetic.index_min());
    double relKinetic = abs((maxKinetic - minKinetic)/maxKinetic);
    REQUIRE(relKinetic < eps);

    vec potential = solver.potentialEnergy.col(0);
    double maxPotential = potential(potential.index_max());
    double minPotential = potential(potential.index_min());
    double relPotential = abs((maxPotential - minPotential)/maxPotential);
    REQUIRE(relPotential < eps);

    mat angular = solver.angularMomentum.col(0);
    double maxAngular = angular(angular.index_max());
    double minAngular = angular(angular.index_min());
    double relAngular = abs((maxAngular - minAngular)/maxAngular);
    REQUIRE(relAngular < eps);
}

TEST_CASE("Elliptic orbit")
//Check that combined kinetic and potential energy is conserved
{
    double eps = 1e-5;

    Planet Earth(vec({1,0,0}), vec({0,M_PI/365.25,0}),1./333000);
    vector<Planet> solarsystem = vector<Planet>{Earth};
    Verlet solver(solarsystem, 1, scale);

    solver.solve(newton, 365.25, 1000000, 1000);

    vec totalEnergy = solver.kineticEnergy.col(0) + solver.potentialEnergy.col(0);
    double maxEnergy = totalEnergy(totalEnergy.index_max());
    double minEnergy = totalEnergy(totalEnergy.index_min());
    double relEnergy = abs((maxEnergy - minEnergy)/maxEnergy);
    REQUIRE(relEnergy < eps);
    cout << relEnergy << endl;
}

TEST_CASE("Three-body")
//Check that combined kinetic and potential energy between both planets
//are conserved
{
    double eps = 1e-5;

    Planet Earth(vec({1,0,0}), vec({0,2*M_PI/365.25,0}),1./333000);
    Planet Jupiter(vec({3,0,0}), vec({0,M_PI/365.25,0}),0.09);
    vector<Planet> solarsystem = vector<Planet>{Earth,Jupiter};
    Verlet solver(solarsystem, 2, scale);

    solver.solve(newton, 365.25, 1000000, 1000);

    vec totalEnergy = solver.kineticEnergy.col(0) + solver.potentialEnergy.col(0)
    + solver.kineticEnergy.col(1) + solver.potentialEnergy.col(1);
    double maxEnergy = totalEnergy(totalEnergy.index_max());
    double minEnergy = totalEnergy(totalEnergy.index_min());
    double relEnergy = abs((maxEnergy - minEnergy)/maxEnergy);
    REQUIRE(relEnergy < eps);

    vec angular = solver.angularMomentum.col(0) + solver.angularMomentum.col(1);
    double maxAngular = angular(angular.index_max());
    double minAngular = angular(angular.index_min());
    double relAngular = abs((maxAngular - minAngular)/maxAngular);
    REQUIRE(relEnergy < eps);
}
