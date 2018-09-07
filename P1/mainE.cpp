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
    Col<double> solution(n);
    Col<double> v(n);
    for(int i=0; i<n; i++)
    {
        v(i) = h*h*100*exp(-10*i*h);
    }

    solution = solve(A,v);
    double* solutionVector = new double[n];
    for(int i = 0; i<n; i++)
    {
        solutionVector[i] = solution(i);
    }

    int N[2] = {1,n};
    FILE * pFile;
    pFile = fopen ("metaArmadillo.bin", "wb");
    fwrite (N, sizeof(int), 2, pFile);
    fclose (pFile);

    //writes data to binary file so python can read it
    pFile = fopen ("myfileArmadillo.bin", "wb");
    fwrite (solutionVector, sizeof(double), n, pFile);
    fclose (pFile);
    return 0;
}
