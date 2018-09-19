#include <iostream>
#include <chrono>
#include <math.h>
#include <armadillo>
#include "func.cpp"
using namespace arma;

int main()
{
    int n = 300;
    double pN = 14;

    double p0 = 0;

    double h = (pN - p0)/n;
    double h2 = h*h;
    Col<double> x = linspace(p0, pN, n);

    mat A = zeros(n, n);
    for(int i = 0; i<n-1; i++)
    {
        A(i,i) = 2/h2+ x(i)*x(i);
        A(i,i+1) = -1/h2;
        A(i+1,i) = -1/h2;
    }
    A(n-1,n-1) = 2/h2 + x(n-1)*x(n-1);

    int k = 0;
    int l = 0;
    double max = 1;

    while (max>1e-10)
    {
        max_element(A, k, l, max, n);
        Jacobi(A, k, l, n);
    }
    //cout << A << endl;

    Col<double> eig = A.diag(0);
    eig = sort(eig);
    cout << eig(0) << endl;
    cout << eig(1) << endl;
    cout << eig(2) << endl;
    cout << eig(3) << endl;
    return 0;
}
