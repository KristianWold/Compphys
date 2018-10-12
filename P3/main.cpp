#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include "solver.h"
#include "planet.h"
using namespace std;
using namespace arma;

double scale = 4*M_PI*M_PI/(365.25*365.25);

vec newton(vec pos, vec vel)
{
    double rCube = pow(norm(pos), 3);
    return -scale/rCube*pos;
}

vec vec3(double x, double y, double z)
{
    vec myVec(3,fill::zeros);
    myVec[0] = x;
    myVec[1] = y;
    myVec[3] = x;
    return myVec;
}

int main(int argc, char const *argv[])
{
    Planet Earth(vec3(1,0,0), vec3(0,2*M_PI/365.25,0),1./333000);
    Planet Jupiter(vec3(3,0,0), vec3(0,M_PI/365.25,0),0.09543);

    vector<Planet> solarsystem = vector<Planet>{Earth,Jupiter};

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

    Verlet solver(solarsystem, 2, scale);
    solver.solve(newton, atof(argv[1]), atof(argv[2]));
    solver.solveEnergy();

    solver.coordinatesToFile("coordinates.txt");
    solver.energyToFile("energy.txt");
    int a = system("python plot.py");

    return 0;
}
