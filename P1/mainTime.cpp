#include <iostream>
#include <math.h>
#include <fstream>
#include "func.cpp"
#include "armadillo"

arma::Col<double> LU(int n, double source(double))
{
    double h = 1/double(n-1);
    arma::Mat<double> A(n,n,arma::fill::zeros);
    A.diag(0).fill(2);
    A.diag(1).fill(-1);
    A.diag(-1).fill(-1);
    arma::Col<double> v(n);
    arma::Col<double> f(n);
    for(int i=0; i<n; i++)
    {
        f(i) = h*h*solution(i*h);
    }
    auto start = std::chrono::high_resolution_clock::now(); //start clock

    v = arma::solve(A,f);

    auto finish = std::chrono::high_resolution_clock::now(); //stop clock
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time for n=" << n <<" using Armadillo LU-decomp: " << elapsed.count() << std::endl;

    return(arma::Col<double>) v;
}

int main(int argc, char const *argv[])
{
    double* v;
    arma::Col<double> V;
    int n = atoi(argv[1]);
    v = solveB(n,source);
    v = solveC(n,source);
    if (n <= 1e5)
    {
        V = LU(n,source);
    }
    return 0;
}
