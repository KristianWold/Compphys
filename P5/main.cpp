#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <armadillo>
#include <fstream>
#include "classes/VMC.hpp"

using namespace std;
using namespace chrono;
using namespace arma;

//Simple trial wave function and accompaying kinetic energy.
//Potential energy is without coloumb term.
//----------------------------------------------------------------------------
inline double acceptAmp1(double* params, double omega, mat &r, vec &delta, int num)
{
    double alpha = params[0];
    double expo = dot(delta, delta+2*r.col(num));
    return exp(-alpha*omega*expo);
}

inline double localKinetic1(double* params, double omega, mat &r)
{
    double R2 = accu(mat(r%r)); // r%r means elementwise multiplication
    double alpha = params[0];
    return -0.5*omega*omega*R2*alpha*alpha + 3*alpha*omega;
}

inline double localPotential1(double* params, double omega, mat &r)
{
    double R2 = accu(mat(r%r)); // r%r means elementwise multiplication
    return 0.5*omega*omega*R2;
}
//----------------------------------------------------------------------------

//Trial wave function with Jastrow factor and accompaying kinetic energy.
//Potential energy is with coloumb term.
//----------------------------------------------------------------------------
inline double acceptAmp2(double* params, double omega, mat &r, vec &delta, int num)
{
    int other;
    if(num == 0) other = 1;
    else other = 0;

    double alpha = params[0];
    double beta = params[1];

    double R1 = -alpha*omega*dot(delta,delta+2*r.col(num));
    double R2 = norm(r.col(num) + delta - r.col(other));
    double R3 = norm(r.col(num) - r.col(other));

    return exp(R1 + R2/(1+beta*R2) - R3/(1+beta*R3));
}

inline double localKinetic2(double* params, double omega, mat &r)
{
    double alpha = params[0];
    double beta = params[1];

    double R12 = norm(r.col(0) - r.col(1));
    double BR = 1 + beta*R12;
    return localKinetic1(params, omega, r) +
    1/(2*BR*BR)*(alpha*omega*R12 - 1/(2*BR*BR) - 2/R12 + 2*beta/BR);
}

inline double localPotential2(double* params, double omega, mat &r)
{
    double R2 = accu(mat(r%r)); // r%r means elementwise multiplication
    return 0.5*omega*omega*R2 + 1/norm(r.col(0) - r.col(1));
}
//----------------------------------------------------------------------------



int main(int argc, char const *argv[]) {

    double* params = new double[1];
    VMC solver1(3, 2, &acceptAmp1, &localKinetic1, &localPotential1);
    ofstream myfile;
    myfile.open("results/noninteracting.txt");
    for(int i=0; i<=20; i++)
    {
        params[0] = 0.5 + 1./20*i;
        Result myResult = solver1.solve(1e5, 1e3, params, 1, false);
        myfile << params[0] << " " << myResult.E << " " << myResult.Var << "\n";
    }

    return 0;
}
