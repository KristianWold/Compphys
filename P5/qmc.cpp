#include<cmath>
#include<random>
#include<iostream>
#include<iomanip>
#include <chrono>

using namespace std;
using namespace chrono;

double acceptAmp(double alpha, double omega, double* delta, double* r)
{
    double expo = delta[0]*(delta[0]+2*r[0]) + delta[1]*(delta[1]+2*r[1]) + delta[2]*(delta[2]+2*r[2]);
    return exp(-alpha*omega*expo);
}

int main(int argc, char const *argv[])
{
    int N = 2;
    int dim = 3;
    double alpha = 0.9;
    double omega = 1;
    int M = 1000;

    double step = 3/sqrt(alpha*omega);
    int accepted = 0;
    double energy = 0;
    double E = 0;
    double E2 = 0;
    double V = 0;
    double R = 0;

    double* delta = new double[dim];
    double** r = new double*[2];
    r[0] = new double[3];
    r[1] = new double[3];

    for (int i=0; i<N; i++)
    {
        for (int j=0; j<dim; j++)
        {
            r[i][j] = 0;
        }
    }

    mt19937_64 engine;
    uniform_real_distribution<double> myRandu = uniform_real_distribution<double>(0,1);

    auto start = high_resolution_clock::now();
    for(int i=0; i<M; i++)
    {
        for(int j=0; j<N; j++)
        {
            delta[0] = (myRandu(engine) - 0.5)*step;
            delta[1] = (myRandu(engine) - 0.5)*step;
            delta[2] = (myRandu(engine) - 0.5)*step;

            if(myRandu(engine) < acceptAmp(alpha, omega, delta, r[j]))
            {
                accepted++;
                r[j][0] += delta[0];
                r[j][1] += delta[1];
                r[j][2] += delta[2];
            }
        }
        R = 0;
        for (int i=0; i<N; i++)
        {
            for (int j=0; j<dim; j++)
            {
                R += r[i][j]*r[i][j];
            }
        }
        energy = 0.5*omega*omega*R*(1-alpha*alpha) + 3*alpha*omega;
        E+= energy;
        E2 += energy*energy;
    }
    auto finish = high_resolution_clock::now();
    E /= M;
    V = E2/M - E*E;

    cout << "Times used: " << duration<double>(finish - start).count() << endl;
    cout << "States accepted: " << accepted << "/" << M << endl;
    cout << "Energy: " << setprecision(5) << E << endl;
    cout << "Variance: " << setprecision(5) << V << endl;
    return 0;
}
