#ifndef FUNC_H
#define FUNC_H

#include <iostream>
#include <math.h>
#include <armadillo>
using namespace arma;
using namespace std;

mat makeMatrix(vec x, double w, vec V(vec x, double w), double h, int n);

void max_element(mat &A, int &k, int &l, double &max, int n);

void jacobi(mat &A, int k, int l, int n);

vec solveJacobi(mat A, int n);

mat solveArma(mat A, int n);

#endif
