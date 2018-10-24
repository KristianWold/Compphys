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
    vec pos_mercury = vec({4.508643080267578E-02, -4.492276414436671E-01, -4.152907939875187E-02}); // [AU]
    vec vel_mercury = vec({2.234601746773343E-02, 4.225668557407336E-03, -1.705402587866491E-03});
    vel_mercury *= 365.242199;
    double M_mercury = 1.65956463e-7; // [solar mass]

    Planet Mercury(pos_mercury, vel_mercury, M_mercury);
    vector<Planet> solarsystem = vector<Planet>{Mercury};
    Solver solver(solarsystem, scale);

    double T = atof(argv[1]);
    int N = atoi(argv[2]);

    solver.solvePerihelion(newton, T, N, "angle1.txt");
    solver.solvePerihelion(einstein, T, N, "angle2.txt");

    system("python3 precession.py");

    return 0;
}
