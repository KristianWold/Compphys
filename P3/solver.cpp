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
               << pos.col(j)(2) << " "

               << vel.col(j)(0) << " "
               << vel.col(j)(1) << " "
               << vel.col(j)(2) << " ";
    }
    myfile << "\n";
}

void Solver::sampleEnergyAndAngular(mat &kinetic, mat &potential, mat &angular, int i)
{
    for(int j=0; j<numPlanets; j++)
    {

        kineticEnergy(i,j) = 0.5*planets[j].M*pow(norm(vel.col(j)),2);

        //potential energy from sun
        potentialEnergy(i,j) = -scale*planets[j].M/norm(pos.col(j));

        //potential energy inbetween planets
        for(int k=j+1; k<numPlanets; k++)
        {
            double temp = -scale*planets[j].M*planets[k].M/
            norm(pos.col(j) - pos.col(k));
            potentialEnergy(i,j) += temp;
            potentialEnergy(i,k) += temp;
        }

        angularMomentum(i,j) = planets[j].M*norm(cross(pos.col(j), vel.col(j)));
    }

}

Solver::Solver(vector<Planet> planets_, double scale_)
{
    planets = planets_;
    numPlanets = planets.size();
    scale = scale_;
}

void Solver::solveVerlet(vec acc(vec, vec), double T, int N, int sampleN)
{
    solved = true;
    double dt = T/N;
    t = linspace(0, T, N);
    pos = zeros(3, numPlanets);
    vel = zeros(3, numPlanets);
    countSample = 0;

    kineticEnergy = zeros(sampleN, numPlanets);
    potentialEnergy = zeros(sampleN, numPlanets);
    angularMomentum = zeros(sampleN, numPlanets);

    for(int i=0; i<numPlanets; i++)
    {
        //initial conditions
        pos.col(i) = planets[i].pos;
        vel.col(i) = planets[i].vel;
    }

    ofstream myfile;
    myfile.open("data.txt");

    mat totalAcc(3, numPlanets, fill::zeros);
    mat prevAcc(3, numPlanets, fill::zeros);
    totalAcceleration(totalAcc, acc);

    for(int i=0; i<N-1; i++)
    {
        if (i%sampleN == 0)
        {
            coordinatesToFile(myfile);
            sampleEnergyAndAngular(kineticEnergy, potentialEnergy,
                angularMomentum, i/sampleN);
        }

        pos = pos + vel*dt + 0.5*totalAcc*dt*dt;
        prevAcc = totalAcc;
        totalAcceleration(totalAcc, acc);
        vel = vel + 0.5*(totalAcc + prevAcc)*dt;


    }
    myfile.close();
}
