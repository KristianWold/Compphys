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

    vec eig;
    int m = 10;
    double average_time = 0;
    for(int n = 10; n <= 100; n += 10)
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

            int k = 0;
            int l = 0;
            double max = 1;

            auto start = high_resolution_clock::now();
            eig = solveJacobi(A, n);
            auto finish = high_resolution_clock::now();

            average_time += duration<double>(finish - start).count();
        }
        average_time = average_time/m;
        myfile << n << " " << average_time << endl;
    }
    myfile.close();
    return 0;
}
