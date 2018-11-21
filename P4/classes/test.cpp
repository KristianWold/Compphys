#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include <stdlib.h>
#include <armadillo>
#include <fstream>
#include <cmath>
#include <map>
#include <random>

#include "metropolis.hpp"
#include "ising.hpp"
using namespace std;

TEST_CASE("Initial energy")
//Check that the energy of an ordered lattice have the correct value
{
    mt19937_64 engine(1111);
    Ising crystal1(ones<Mat<int> >(2,2), 2, 1, 1, engine);
    REQUIRE(crystal1.energy == -8);

    Ising crystal2(ones<Mat<int> >(4,4), 4, 1, 1, engine);
    REQUIRE(crystal2.energy == -32);
}

TEST_CASE("Periodic boundary")
//Check that the energy of an ordered lattice have the correct value
{
    mt19937_64 engine(1111);
    Ising crystal2(ones<Mat<int> >(4,4), 4, 1, 1, engine);
    REQUIRE(crystal2.periodic(4) == 0);
    REQUIRE(crystal2.periodic(-1) == 3);
}
