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


void Solver::totalAcceleration(mat &totalAcc, vec acc(vec, vec))
{
    totalAcc = zeros(3, numPlanets);
    for(int j=0; j<numPlanets; j++)
    {
        totalAcc.col(j) += acc(pos.col(j), vel.col(j));
        for(int k=j+1; k<numPlanets; k++)
        {
            vec relpos = pos.col(j) - pos.col(k);
            vec relvel = vel.col(j) - vel.col(k);
            mat temp = acc(relpos, relvel);
            totalAcc.col(j) += planets[k].M*temp;
            totalAcc.col(k) -= planets[j].M*temp;
        }
    }
}

void Solver::coordinatesToFile(ofstream &myfile)
{
    for(int j=0; j<numPlanets; j++)
    {
        myfile << pos.col(j)(0) << " "
               << pos.col(j)(1) << " "
               << pos.col(j)(2) << " ";
    }
    myfile << "\n";
}

void Solver::sampleEnergyAndAngular(mat &kinetic, mat &potential,
                                    vec &energyAllPlanets, vec &angular, int i)
{
    vec totalAngular(3,fill::zeros);
    for(int j=0; j<numPlanets; j++)
    {

        kineticEnergy(i,j) = 0.5*planets[j].M*pow(norm(vel.col(j)),2);

        //potential energy from sun
        potentialEnergy(i,j) = -scale*planets[j].M/norm(pos.col(j));
        energyAllPlanets(i) = +potentialEnergy(i,j);
        //potential energy inbetween planets
        for(int k=j+1; k<numPlanets; k++)
        {
            double temp = -scale*planets[j].M*planets[k].M/
            norm(pos.col(j) - pos.col(k));
            potentialEnergy(i,j) += temp;
            potentialEnergy(i,k) += temp;
            energyAllPlanets(i)  += temp;
        }
        energyAllPlanets(i) += kineticEnergy(i,j);
        totalAngular += planets[j].M*cross(pos.col(j), vel.col(j));
    }
    angularMomentum(i) = norm(totalAngular);

}

void Solver::euler(vec acc(vec, vec))
{
    prevVel = vel;
    vel = vel + totalAcc*dt;
    pos = pos + prevVel*dt;
    totalAcceleration(totalAcc, acc);
}

void Solver::verlet(vec acc(vec, vec))
{
    pos = pos + vel*dt + 0.5*totalAcc*dt*dt;
    prevAcc = totalAcc;
    totalAcceleration(totalAcc, acc);
    vel = vel + 0.5*(totalAcc + prevAcc)*dt;
}

Solver::Solver(vector<Planet> planets_, double scale_)
{
    planets = planets_;
    numPlanets = planets.size();
    scale = scale_;
}

void Solver::solve(int method, vec acc(vec, vec), double T, int N, int sampleN)
{
    solved = true;
    dt = T/N;
    //t = linspace(0, T, N);
    pos = zeros(3, numPlanets);
    vel = zeros(3, numPlanets);

    kineticEnergy      = zeros(N/sampleN, numPlanets);
    potentialEnergy    = zeros(N/sampleN, numPlanets);
    energyAllPlanets   = zeros(N/sampleN);
    angularMomentum    = zeros(N/sampleN);

    for(int i=0; i<numPlanets; i++)
    {
        //initial conditions
        pos.col(i) = planets[i].pos;
        vel.col(i) = planets[i].vel;
    }

    ofstream myfile;
    myfile.open("data.txt");

    totalAcc = zeros(3, numPlanets);
    prevAcc = zeros(3, numPlanets);
    prevVel = zeros(3, numPlanets);

    totalAcceleration(totalAcc, acc);

    coordinatesToFile(myfile);
    sampleEnergyAndAngular(kineticEnergy, potentialEnergy,
        energyAllPlanets, angularMomentum, 0);

    for(int i=0; i<N-1; i++)
    {
        if(method == 1) euler(acc);
        if(method == 2) verlet(acc);
        if ((i+1)%sampleN == 0)
        {
            coordinatesToFile(myfile);
            sampleEnergyAndAngular(kineticEnergy, potentialEnergy,
                energyAllPlanets, angularMomentum, (i+1)/sampleN);
        }

    }
    myfile.close();
}

double Solver::kineticFluctuation(int i)
{
    vec kinetic = kineticEnergy.col(i);
    double maxKinetic = kinetic(kinetic.index_max());
    double minKinetic = kinetic(kinetic.index_min());
    return abs((maxKinetic - minKinetic)/maxKinetic);
}

double Solver::potentialFluctuation(int i)
{
    vec potential = potentialEnergy.col(i);
    double maxPotential = potential(potential.index_max());
    double minPotential = potential(potential.index_min());
    return abs((maxPotential - minPotential)/maxPotential);
}

double Solver::totalEnergyFluctuation()
{

    double maxEnergy = energyAllPlanets(energyAllPlanets.index_max());
    double minEnergy = energyAllPlanets(energyAllPlanets.index_min());
    return abs((maxEnergy - minEnergy)/maxEnergy);
}

double Solver::angularFluctuation()
{
    double maxAngular = angularMomentum(angularMomentum.index_max());
    double minAngular = angularMomentum(angularMomentum.index_min());
    return abs((maxAngular - minAngular)/maxAngular);
}
