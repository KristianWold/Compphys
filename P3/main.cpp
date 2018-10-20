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
    Planet Earth(vec({1,0,0}), vec({0,2*M_PI,0}),1./333000);
    vector<Planet> solarsystem = vector<Planet>{Earth};
    Solver solverEuler(solarsystem, scale);
    Solver solverVerlet(solarsystem, scale);

    double T = atof(argv[1]);
    int N = atoi(argv[2]);
    int sampleN = atoi(argv[3]);

    solverEuler.solve(1, newton, T, N, sampleN);
    system("python3 plot.py pos");

    solverVerlet.solve(2, newton, T, N, sampleN);

    system("python3 plot.py pos");

    ofstream myfile;
    myfile.open("energy.txt");
    for(int i=0; i<N/sampleN; i++)
    {
        myfile << setprecision(8)
               << solverEuler.energyAllPlanets(i) << " "
               << solverVerlet.energyAllPlanets(i) << "\n";
    }
    myfile.close();
    system("python3 plot.py energy");

    myfile.open("fluctuation.txt");

    for(int n = 100; n<=1e8; n*=10)
    {
        Solver solver(solarsystem, scale);
        solver.solve(1, newton, T, n, n/10);
        myfile << n << " " << solver.totalEnergyFluctuation() << "\n";
        cout << n << endl;
    }
    myfile.close();

    system("python3 plot.py fluctuation");

    myfile.open("fluctuation.txt");

    for(int n = 100; n<=1e8; n*=10)
    {
        Solver solver(solarsystem, scale);
        solver.solve(2, newton, T, n, n/10);
        myfile << n << " " << solver.totalEnergyFluctuation() << "\n";
        cout << n << endl;
    }
    myfile.close();

    system("python3 plot.py fluctuation");

    return 0;
}



    /*int count = 0;
    ifstream myfile;
    myfile.open("init.txt");
    string line;
    while(getline(myfile, line))
    {

        istringstream buf(line);
        istream_iterator<string> beg(buf), end;
        vector<string> tokens(beg, end);

        solarsystem[count].pos(0) = stof(tokens[0]);
        solarsystem[count].pos(1) = stof(tokens[1]);
        solarsystem[count].pos(2) = stof(tokens[2]);

        solarsystem[count].vel(0) = stof(tokens[3]);
        solarsystem[count].vel(1) = stof(tokens[4]);
        solarsystem[count].vel(2) = stof(tokens[5]);

        count++;
    }
    myfile.close();
    */
