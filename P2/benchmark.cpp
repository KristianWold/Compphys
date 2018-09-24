#include <iostream>
#include <armadillo>
#include <chrono>
#include <fstream>

#include "func.h"

using namespace arma;
using namespace std;
using namespace chrono;

int main()
{
    ofstream myfile;
    myfile.open("benchmark.txt");

    vec eigval;
    mat eigvec;
    int m = 10;
    double average_time_jacobi = 0;
    double average_time_arma = 0;
    for(int n = 10; n <= 120; n += 10)
    {

        for(int i = 0; i < m; i++)
        {
            mat A = zeros(n, n);
            for (int i = 0; i < n-1; i++)
            {
                A(i,i) = 2;
                A(i,i+1) = -1;
                A(i+1,i) = -1;
            }
            A(n-1,n-1) = 2;

            auto start = high_resolution_clock::now();
            eigval = solveJacobi(A, n);
            auto finish = high_resolution_clock::now();
            average_time_jacobi += duration<double>(finish - start).count();

            start = high_resolution_clock::now();
            eigvec = solveArma(A, n);
            finish = high_resolution_clock::now();
            average_time_arma += duration<double>(finish - start).count();
        }
        average_time_jacobi = average_time_jacobi/m;
        average_time_arma = average_time_arma/m;
        myfile << n << " " << average_time_jacobi << " " << average_time_arma << endl;
    }
    myfile.close();
    return 0;
}
