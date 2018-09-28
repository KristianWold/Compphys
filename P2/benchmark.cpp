#include <iostream>
#include <armadillo>
#include <chrono>
#include <fstream>

#include "jacobi.h"

using namespace arma;
using namespace std;
using namespace chrono;

int main()
{
    ofstream myfile;
    myfile.open("benchmark.txt");

    vec eigval;
    mat eigvec;
    int N;
    int m = 1;
    double average_time_jacobi = 0;
    double average_time_arma = 0;
    //Runs Jacobis method and Armadillos eig_sym n = 10, 20, ..., 120
    for(int n = 10; n <= 400; n += 20)
    {
        //Runs the methods m = 10 times for each n, then averages the times
        for(int i = 0; i < m; i++)
        {
            mat A = zeros(n, n);
            for (int j = 0; j < n-1; j++)
            {
                A(j,j) = 2;
                A(j,j+1) = -1;
                A(j+1,j) = -1;
            }
            A(n-1,n-1) = 2;

            auto start = high_resolution_clock::now();
            N = solveJacobi(A, eigval, eigvec, n);
            auto finish = high_resolution_clock::now();
            average_time_jacobi += duration<double>(finish - start).count();

            start = high_resolution_clock::now();
            eig_sym(eigval, eigvec, A);
            finish = high_resolution_clock::now();
            average_time_arma += duration<double>(finish - start).count();
        }
        //Averages the times
        average_time_jacobi = average_time_jacobi/m;
        average_time_arma = average_time_arma/m;
        //Writes to file
        myfile << n << " " << average_time_jacobi << " " << average_time_arma
        <<" " << N << endl;
    }
    myfile.close();
    return 0;
}
