#ifndef FUNC_H
#define FUNC_H

#include <iostream>
#include <math.h>
#include <armadillo>
using namespace arma;
using namespace std;

void max_element(mat, int, int, double, int);

void Jacobi(mat, int, int, int)

vec SolveJacobi(double, vec, vec, double, int)

mat SolveArma(int, double, vec, vec, double w)

#endif
