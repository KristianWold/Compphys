#include <iostream>
#include <iomanip>
#include <armadillo>
#include "func.cpp"
#include <math.h>
using namespace arma;
using namespace std;

int main()
{

    //test max_element
    bool passed = true;

    int n = 4;
    mat A = zeros(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A(i,j) = 10*i + j;
        }
    }

    int k = 0;
    int l = 0;
    double max = 0;

    max_element(A, k, l, max, n);
    if (k != 2 or l != 3){passed = false;}

    A(1,2) = -100;
    max = 0;

    max_element(A, k, l, max, n);
    if (k != 1 or l != 2){passed = false;}

    if (passed){cout << "Found larget element" << endl;}
    else{cout << "Didn't find larget element" << endl;}

    //test eigenvalues
    passed = true;

    n = 10;
    A = zeros(n, n);
    for (int i = 0; i < n-1; i++)
    {
        A(i,i) = 2;
        A(i,i+1) = -1;
        //A(i+1,i) = -1;
    }
    A(n-1,n-1) = 2;

    k = 0;
    l = 0;
    max = 1;
    while (max>1e-10)
    {
        max_element(A, k, l, max, n);
        Jacobi(A, k, l, n);
    }

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
    if (passed){cout << "Correct eigenvalues" << endl;}
    else{cout << "Incorrect eigenvalues" << endl;}
    return 0;
}
