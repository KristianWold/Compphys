#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include "Classes/solver.h"
#include "Classes/planet.h"
using namespace std;
using namespace arma;

double scale = 4*M_PI*M_PI;

vec newton(vec pos, vec vel)
{
    double rCube = pow(norm(pos), 3);
    return -scale/rCube*pos;
}

int main(int argc, char const *argv[])
{
    Planet Earth(vec({1,0,0}), vec({0,2*M_PI,0}),3e-6);
    Planet Jupiter(vec({5,0,0}), vec({0,3,0}), 0.001);
    vector<Planet> solarsystem = vector<Planet>{Earth, Jupiter};

    double T = 10;
    int N = 1000000;
    int sampleN = 10;

    vector<double> mass{0.001, 0.01, 1};

    ofstream myfile;
    myfile.open("fluctuation.txt");

    for(int n = 100; n<=1e8; n*=10)
    {
        Solver solver(solarsystem, scale);
        solver.solve(1, newton, 5, n, n/100, "data.txt");
        myfile << n << " " << solver.totalEnergyFluctuation() << "\n";
        cout << n << endl;
    }
    myfile.close();

    system("python3 plot.py fluctuation EarthAndJupiter");
    myfile.open("fluctuation.txt");

    for(int i = 1; i<=3; i++)
    {
        solarsystem[1].M = mass[i-1];
        Solver solver(solarsystem, scale);
        solver.solve(2, newton, T, N, sampleN, "data" + to_string(i) + ".txt");
    }

    system("python3 plot.py earthAndJupiter");



    return 0;
}
