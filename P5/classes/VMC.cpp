#include<cmath>
#include<random>
#include<iostream>
#include<iomanip>
#include<chrono>
#include<armadillo>
#include<fstream>

using namespace std;
using namespace chrono;
using namespace arma;

//function pointer type to pass special functons to the class VMC
typedef double (*myfunc1)(double*, double, mat&, vec&, int);
typedef double (*myfunc2)(double*, double, mat&);

struct result
{
    double E;
    double Var;
    double kinetic;
    double potential;
    double R12;
    int accepted;
};

class VMC
{
public:
    myfunc1 acceptAmp;
    myfunc2 localKinetic, localPotential;
    int numDim, numParticles, numCycles, accepted;
    double kineticE, potentialE, k_E, p_E, E, E2, Var, R12, step;
    mat pos;
    vec delta;

    mt19937_64 engine;
    uniform_real_distribution<double> myRandu = uniform_real_distribution<double>(0,1);

    VMC(int numDim, int numParticles, myfunc1 acceptAmp, myfunc2 localKinetic, myfunc2 localPotential)
    {
        this->numDim = numDim;
        this->numParticles = numParticles;
        this->acceptAmp = acceptAmp;
        this->localKinetic = localKinetic;
        this->localPotential = localPotential;

        pos = zeros(numDim, numParticles);
        delta = zeros(numDim);
    }

    result solve(int numCycles, int preCycles, double* params, double omega)
    {
        ofstream myfile;
        step = 2.1/sqrt(params[0]*omega);
        kineticE = 0;
        potentialE = 0;
        k_E = 0;
        p_E = 0;
        E = 0;
        E2 = 0;
        Var = 0;
        R12 = 0;
        accepted = 0;

        myfile.open("data.txt");
        for(int i=0; i<numCycles+preCycles; i++)
        {
            for(int j=0; j<numParticles; j++)
            {
                for(int k=0; k<numDim; k++)
                {
                    delta(k) = (myRandu(engine) - 0.5)*step;
                }
                if(myRandu(engine)<acceptAmp(params, omega, pos, delta, j))
                {
                    pos.col(j) += delta;
                    accepted++;
                }
            }
            if(i>=preCycles)
            {
                kineticE = localKinetic(params, omega, pos);
                potentialE = localPotential(params, omega, pos);
                myfile << pos(0,0) << " " << pos(1,0) << " " << pos(2,0) << " "
                       << pos(0,1) << " " << pos(1,1) << " " << pos(2,1) << " "
                       << kineticE + potentialE << "\n";

                k_E += kineticE;
                p_E += potentialE;
                E += kineticE + potentialE;
                E2 += (kineticE + potentialE)*(kineticE + potentialE);
                R12 += norm(pos.col(0) - pos.col(1));
            }
        }
        myfile.close();

        E /= numCycles;
        E2 /= numCycles;
        Var = E2 - E*E;
        k_E /= numCycles;
        p_E /= numCycles;
        R12 /= numCycles;

        result myResult = {E, Var, k_E, p_E, R12, accepted};
        return myResult;
    }
};

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
    int numCycles = atoi(argv[1]);
    int preCycles = atoi(argv[2]);
    double omega = atof(argv[3]);
    double alpha = atof(argv[4]);
    double beta = atof(argv[5]);

    double* params = new double[2];
    params[0] = alpha;
    params[1] = beta;

    VMC solver1(3, 2, &acceptAmp2, &localKinetic2, &localPotential2);
    result myResult = solver1.solve(numCycles, preCycles, params, omega);

    cout << "Energy = " << myResult.E << endl;
    cout << "Variance = " << myResult.Var << endl;
    cout << "Ratio = " << myResult.kinetic/myResult.potential << endl;
    cout << "<|r1 - r2|> = " << myResult.R12 << endl;
    cout << "Accepted " << myResult.accepted << endl;


    return 0;
}
