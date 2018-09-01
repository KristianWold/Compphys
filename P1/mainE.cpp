#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;

int main(int argc, char const *argv[])
{
    int n = atof(argv[1]);
    double h = 1./(n+1);
    Mat<double> A(n,n);
    A.zeros();
    for(int i=0; i<n-1; i++)
    {
        A(i,i) = 2.;
        A(i+1,i) = -1.;
        A(i,i+1) = -1.;
    }
    A(n-1,n-1) = 2.;
    Col<double> x(n);
    Col<double> v(n);
    for(int i=0; i<n-1; i++)
    {
        v(i) = 100*exp(-10*i*h);
    }

    x = solve(A, v);


    return 0;
}
