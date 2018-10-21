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

    solverEuler.solve(1, newton, 5, 100, 1);
    system("python3 plot.py singlePlanet Euler 5 100");

    solverEuler.solve(1, newton, 5, 1000, 1);
    system("python3 plot.py singlePlanet Euler 5 1000");

    solverEuler.solve(1, newton, 5, 10000, 1);
    system("python3 plot.py singlePlanet Euler 5 10000");


    solverVerlet.solve(2, newton, 5, 100, 1);
    system("python3 plot.py singlePlanet Verlet 5 100");

    solverVerlet.solve(2, newton, 5, 1000, 1);
    system("python3 plot.py singlePlanet Verlet 5 1000");

    solverVerlet.solve(2, newton, 5, 10000, 1);
    system("python3 plot.py singlePlanet Verlet 5 10000");


    ofstream myfile;
    myfile.open("energy.txt");
    for(int i=0; i<10000; i++)
    {
        myfile << setprecision(8)
               << solverEuler.energyAllPlanets(i) << " "
               << solverVerlet.energyAllPlanets(i) << "\n";
    }
    myfile.close();
    system("python3 plot.py energy 5 10000");

    myfile.open("fluctuation.txt");

    for(int n = 100; n<=1e8; n*=10)
    {
        Solver solver(solarsystem, scale);
        solver.solve(1, newton, 5, n, n/100);
        myfile << n << " " << solver.totalEnergyFluctuation() << "\n";
        cout << n << endl;
    }
    myfile.close();

    system("python3 plot.py fluctuation Euler");

    myfile.open("fluctuation.txt");

    for(int n = 100; n<=1e8; n*=10)
    {
        Solver solver(solarsystem, scale);
        solver.solve(2, newton, 5, n, n/10);
        myfile << n << " " << solver.totalEnergyFluctuation() << "\n";
        cout << n << endl;
    }
    myfile.close();

    system("python3 plot.py fluctuation Verlet");

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
