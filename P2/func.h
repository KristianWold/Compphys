#ifndef FUNC_H
#define FUNC_H

#include <iostream>
#include <math.h>
#include <armadillo>
using namespace arma;
using namespace std;

mat makeMatrix(vec V, double h, int n);

void max_element(mat &A, int &k, int &l, double &max, int n);

void jacobi(mat &A, mat &R, int k, int l, int n);

int solveJacobi(mat A, vec &eigval, mat &eigvec, int n);

#endif
