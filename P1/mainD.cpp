#include <iostream>
#include <math.h>
#include <fstream>
#include "func.cpp"
using namespace std;

int main()
{
    int N = 30;
    double *err = vector(N);
    double *logh = vector(N);
    //checking error for n=10 to n=10^7
    for(int i=0; i<=N; i++)
    {
        int n = int(pow(10, 1 + 6*i/30.));
        double h = 1./(n + 1);
        double *v = solveC(n);
        double maxerr = 0.;
        for(int j=0; j<n; j++)
        {
            double u = solution(j*h);
            //calculating the error
            double e = abs((v[j] - u)/u);
            //finding the maximum error
            if (e > maxerr)
            {
                maxerr = e;
            }
        }
        err[i] = maxerr;
        logh[i] = h;
        free_vector(v);
    }
    FILE * pFile;
    pFile = fopen ("error.bin", "wb");
    fwrite (logh, sizeof(double), N, pFile);
    fwrite (err, sizeof(double), N, pFile);
    fclose (pFile);

    return 0;
}
