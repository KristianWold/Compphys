#include <iostream>
#include <chrono>
#include <math.h>
#include <armadillo>
#include "func.cpp"
using namespace arma;

int main()
{
    int n = 5;
    mat A = zeros(n, n);
    for (int i = 0; i < n-1; i++)
    {
        A(i,i) = 2;
        A(i,i+1) = -1;
        A(i+1,i) = -1;
    }
    A(n-1,n-1) = 2;

    int k = 0;
    int l = 0;
    double max = 1;
    while (max>1e-10)
    {
        max_element(A, k, l, max, n);
        Jacobi(A, k, l, n);
    }
    cout << A << endl;
    return 0;
}
