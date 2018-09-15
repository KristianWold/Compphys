#include <iostream>
#include <chrono>
#include <math.h>
#include <armadillo>
using namespace arma;

void Jacobi(mat<double> A, int k, int l)
{
    double tau = (A(l)(l) - A(k)(k))/(2*A(k)(l));
    double t = 1/(tau + (1 - 2*signbit(tau))*sqrt(1 + tau*tau));
    double c = 1/sqrt(1 + t*t);
    double s = t*c;

    for (int i = 0; i<l; i++)
    {
        A(i)(k) = A(i)(k)*c - A(i)(l)*s;
    }

    for (int i = k+1; i<n; i++)
    {
        A(k)(i) = A(i)(k)*c - A(l)(i)*s;
    }

    for (int i = 0; i<n-; i++)
    {
        A(i)(k) = A(i)(k)*c - A(i)(l)*s;
    }

    for (int i = k+1; i<n; i++)
    {
        A(k)(i) = A(i)(k)*c - A(l)(i)*s;
    }

    for (int i = b+1; i<a; i++)
    {
        double tempele = B(k)(i) = B(i)(k)*c - B(i)(l)*s;
        B(i)(l) = B(l)(i) = B(i)(l)*c + B(i)(k)*s;
        B(i)(k) = tempele;
    }

    for (int i = a+1; i<n; i++)
    {
        double tempele = B(k)(i) = B(i)(k)*c - B(i)(l)*s;
        B(i)(l) = B(l)(i) = B(i)(l)*c + B(i)(k)*s;
        B(i)(k) = tempele;
    }

    double cs = c*s;
    double Acs = 2*B(k)(l)*cs;
    double c2 = pow(c,2);
    double s2 = pow(s,2);

    double tempele1 = B(k)(k)*c2 - Acs + B(l)(l)*s2;
    double tempele2 = B(l)(l)*c2 + Acs + B(k)(k)*s2;
    B(k)(l) = B(l)(k) = (B(k)(k) - B(l)(l))*cs + B(k)(l)*(c2-s2);
    B(k)(k) = tempele1;
    B(l)(l) = tempele2;
}
