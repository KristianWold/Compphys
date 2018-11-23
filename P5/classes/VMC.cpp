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

class VMC
{
public:
    myfunc1 acceptAmp;
    myfunc2 localEnergy;
    int numDim, numParticles, numCycles, accepted;
    double energy, E, E2, V, step;
    mat pos;
    vec delta;

    mt19937_64 engine;
    uniform_real_distribution<double> myRandu = uniform_real_distribution<double>(0,1);

    VMC(int numDim, int numParticles, myfunc1 acceptAmp, myfunc2 localEnergy)
    {
        this->numDim = numDim;
        this->numParticles = numParticles;
        this->acceptAmp = acceptAmp;
        this->localEnergy = localEnergy;

        pos = zeros(numDim, numParticles);
        delta = zeros(numDim);
    }

    void solve(int numCycles, int preCycles, double* params, double omega)
    {
        ofstream myfile;
        step = 2.3/sqrt(params[0]*omega);
        energy = 0;
        E = 0;
        E2 = 0;
        V = 0;
        accepted = 0;
        myfile.open("data.txt");
        auto start = high_resolution_clock::now();
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
                    accepted ++;
                }
            }
            if(i>=preCycles)
            {
                myfile << pos(0,0) << " " << pos(1,0) << " " << pos(2,0) << " "
                       << pos(0,1) << " " << pos(1,1) << " " << pos(2,1) << "\n";
                energy = localEnergy(params, omega, pos);
                E += energy;
                E2 += energy*energy;
            }
        }
        auto finish = high_resolution_clock::now();
        myfile.close();

        E /= numCycles;
        E2 /= numCycles;
        V = E2 - E*E;
        //cout << "Times used: " << duration<double>(finish - start).count() << endl;
        cout << "Energy: " << E << "  For alpha=" << params[0] << " ,beta=" << params[1] << endl;
        cout << "Variance: " << V << endl;
        //cout << "Accepted moves: " << accepted << "/" << (numParticles*numCycles) << endl;

    }
};


inline double acceptAmp1(double* params, double omega, mat &r, vec &delta, int num)
{
    double alpha = params[0];
    double expo = dot(delta, delta+2*r.col(num));
    return exp(-alpha*omega*expo);
}

inline double localEnergy1(double* params, double omega, mat &r)
{
    double R2 = accu(mat(r%r));
    double alpha = params[0];
    return 0.5*omega*omega*R2*(1-alpha*alpha) + 3*alpha*omega + 1/norm(r.col(0) - r.col(1));
}


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

inline double localEnergy2(double* params, double omega, mat &r)
{
    double alpha = params[0];
    double beta = params[1];

    double R12 = norm(r.col(0) - r.col(1));
    double BR = 1 + beta*R12;
    return localEnergy1(params, omega, r) +
    1/(2*BR*BR)*(alpha*omega*R12 - 1/(2*BR*BR) - 2/R12 + 2*beta/BR);
}

int main(int argc, char const *argv[]) {
    int numCycles = atoi(argv[1]);
    int preCycles = atoi(argv[2]);
    double omega = atof(argv[3]);
    double alpha = atof(argv[4]);
    double beta = atof(argv[5]);

    double* params = new double[2];
    params[0] = alpha;
    params[1] = beta;
    VMC solver(3, 2, &acceptAmp2, &localEnergy2);
    solver.solve(numCycles, preCycles, params, omega);
    /*
    for(int i=0; i<50; i++)
    {
        for(int j=0; j<50; j++)
        {
            params[0] =  0.75+0.01*i;
            params[1] =  0.75+0.01*j;

            solver.solve(numCycles, preCycles, params, omega);

        }
    }*/
    return 0;
}
