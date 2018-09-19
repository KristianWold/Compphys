#include <iostream>
#include <chrono>
#include <math.h>
#include <armadillo>
using namespace arma;

void max_element(mat &A, int &k, int &l, double &max, int n)
{
    max = 0;
    for(int i = 0; i<n; i++)
    {
        for(int j = i+1; j<n; j++)
        {
            if (abs(A(i,j)) > max)
            {
                max = abs(A(i,j));
                k = i;
                l = j;
            }
        }
    }
}

void Jacobi(mat &A, int k, int l, int n)
{
    double Akk = A(k,k);
    double All = A(l,l);
    double Akl = A(k,l);

    double t;
    double tau = (All - Akk)/(2*Akl);
    if (tau>=0){t = 1/(tau + sqrt(1 + tau*tau));}
    else {t = -1/(tau - sqrt(1 + tau*tau));}

    double c = 1/sqrt(1 + t*t);
    double s = t*c;

    for (int i = 0; i<n; i++)
    {
        double temp = A(i,k);
        A(i,k) = A(i,k)*c - A(i,l)*s;
        A(i,l) = A(i,l)*c + temp*s;
    }

    for (int i = 0; i<n; i++)
    {
        double temp = A(k,i);
        A(k,i) = A(k,i)*c - A(l,i)*s;
        A(l,i) = A(l,i)*c + temp*s;
    }

    double c2 = c*c;
    double s2 = s*s;

    A(k,k) = Akk*c2 - 2*Akl*c*s + All*s2;
    A(l,l) = All*c2 + 2*Akl*c*s + Akk*s2;
    A(k,l) = 0;
}
