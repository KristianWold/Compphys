#include <iostream>
#include <math.h>
#include <fstream>
#include "func.cpp"
using namespace std;

int main()
{
    int N = 30;
    double *maxerror = vector(N);
    double *h = vector(N);
    //checking error for n=10 to n=10^7
    for(int i=0; i<N; i++)
    {
        int n = int(pow(10, 1 + 6*i/(N-1.)));
        double *v = solveC(n);
        h[i] = 1/(float(n) + 1);
        double u = solution(n/2*h[i]);

        maxerror[i] = abs((v[n/2] - u)/u);
        /*
        for(int j=1; j<n-1; j++)
        {
            double u = solution(j*h[i]);
            //calculating the error
            double error = abs((v[j] - u)/u);
            //finding the maximum error
            if (error > maxerror[i])
            {
                maxerror[i] = error;
            }
        }
        */
        free_vector(v);
    }
    FILE * pFile;
    pFile = fopen ("error.bin", "wb");
    fwrite (h, sizeof(double), N, pFile);
    fwrite (maxerror, sizeof(double), N, pFile);
    fclose (pFile);

    return 0;
}
