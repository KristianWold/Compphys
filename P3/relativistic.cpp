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
double c = 63239.7263;          //AU/yr

inline vec newton(vec pos, vec vel)
{
    double r = norm(pos);
    return -scale/pow(r, 3)*pos;
}

inline vec einstein(vec pos, vec vel)
{
    double r = norm(pos);
    double l = abs(norm(cross(pos, vel)));
    return -scale/pow(r, 3)*(1 + 3*l*l/(r*r*c*c))*pos;
}

int main(int argc, char const *argv[])
{
    vec pos_mercury = vec({-1.323E-01, -4.393E-01, -2.444E-02}); // [AU]
    vec vel_mercury = vec({2.132E-02, 6.577E-03, -2.494E-03});
    vel_mercury *= 365.25;
    double M_mercury = 0.166E-6; // [solar mass]

    Planet Mercury(pos_mercury, vel_mercury, M_mercury);
    vector<Planet> solarsystem = vector<Planet>{Mercury};
    Solver solver(solarsystem, scale);

    double T = atof(argv[1]);
    int N = atoi(argv[2]);
    int sampleN = atoi(argv[3]);

    solver.solve(2, einstein, T, N, sampleN);

    system("python3 precession.py");

    return 0;
}
