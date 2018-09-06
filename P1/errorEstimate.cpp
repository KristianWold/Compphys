#include <iostream>
#include <math.h>
#include <fstream>
#include "func.cpp"
using namespace std;

int main()
{
    int N = 7;
    double *h = vector(N);
    double *maxerror = vector(N);

    //checking error for n=10 to n=10^7
    for(int i=0; i<N; i++)
    {
        int n = pow(10, i+1);
        //Numerical solution
        double *v = solveC(n,source);
        h[i] = 1/(double(n)-1);

        maxerror[i] = 0;
        //Truncates the function on each side to avoid division by zero
        for(int j=n/10; j<9*n/10; j++)
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
        cout << maxerror[i] << endl;

        free_vector(v);
    }
    FILE * pFile;
    pFile = fopen ("error.bin", "wb");
    fwrite(h, sizeof(double), N, pFile);
    fwrite(maxerror, sizeof(double), N, pFile);
    fclose(pFile);
    return 0;
}
