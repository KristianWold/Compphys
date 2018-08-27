#include <iostream>
#include <math.h>
#include <fstream>
#include "func.cpp"
using namespace std;

int main()
{
    int N = 7;
    //checking error for n=10, 10^2, ..., 10^N
    for(int i=1; i<=N; i++)
    {
        int n = pow(10,i);
        double h = 1./(n + 1);
        double *v = solveC(n);
        double maxerr = 0.;
        for(int j=0; j<n; j++)
        {
            double u = solution(j*h);
            //calculating the error
            double e = log10(abs((v[j] - u)/u));
            //finding the maximum error
            if (e > maxerr)
            {
                maxerr = e;
            }
        }
        cout<< maxerr << endl;
        free_vector(v);
    }

    return 0;
}
