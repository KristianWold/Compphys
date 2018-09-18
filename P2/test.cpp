#include <iostream>
#include <iomanip>
#include <armadillo>
#include "func.cpp"
#include <math.h>
using namespace arma;
using namespace std;

int main()
{
    int n = 100;
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

    bool passed = true;
    double eps = 1e-6;
    double analytical;
    Col<double> eig = A.diag(0);
    eig = sort(eig);
    for (int i = 0; i < n; i++)
    {
        analytical = 2*(1 - cos((i+1)*M_PI/(n+1)));
        if (abs(eig(i) - analytical)> eps)
        {
          passed = false;
        }
        //cout << "Numerical: "  << setprecision(3) << eig(i);
        //cout << "  Analytical: " << setprecision(3) << analytical << endl;
    }
    if (passed){cout << "Test passed" << endl;}
    else{cout << "Test failed" << endl;}
    return 0;
}
