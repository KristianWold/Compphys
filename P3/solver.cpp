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

void Verlet::totalAcceleration(mat &totalAcc, vec acc(vec, vec), int i)
{
    totalAcc = zeros(3, numPlanets);
    for(int j=0; j<numPlanets; j++)
    {
        totalAcc.col(j) += acc(pos.slice(i).col(j), vel.slice(i).col(j));
        for(int k=j+1; k<numPlanets; k++)
        {
            vec relpos = pos.slice(i).col(j) - pos.slice(i).col(k);
            vec relvel = vel.slice(i).col(j) - vel.slice(i).col(k);
            mat temp = acc(relpos, relvel);
            totalAcc.col(j) += planets[k].M*temp;
            totalAcc.col(k) -= planets[j].M*temp;
        }
    }
}

Verlet::Verlet(vector<Planet> p, int n, double s)
{
    planets = p;
    numPlanets = n;
    scale = s;
}

void Verlet::solve(vec acc(vec, vec), double T, double dt)
    {
        solved = true;
        N = int(T/dt);
        t = linspace(0, T, N);
        pos = zeros(3, numPlanets, N);
        vel = zeros(3, numPlanets, N);
        for(int i=0; i<numPlanets; i++)
        {
            pos.slice(0).col(i) = planets[i].pos;
            vel.slice(0).col(i) = planets[i].vel;
        }

        mat totalAcc(3, numPlanets, fill::zeros);
        mat prevAcc(3, numPlanets, fill::zeros);
        totalAcceleration(totalAcc, acc, 0);
        for(int i=0; i<N-1; i++)
        {
            pos.slice(i+1) = pos.slice(i) + vel.slice(i)*dt + 0.5*totalAcc*dt*dt;
            prevAcc = totalAcc;
            totalAcceleration(totalAcc, acc, i+1);
            vel.slice(i+1) = vel.slice(i) + 0.5*(totalAcc + prevAcc)*dt;
        }
}

void Verlet::solveEnergy()
{
    if (solved == false)
    {
        throw invalid_argument("Must run .solve() first");
    }
    solvedEnergy = true;
    kineticEnergy = zeros(N, numPlanets);
    potentialEnergy = zeros(N, numPlanets);
    angularMomentum = zeros(N, numPlanets);

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<numPlanets; j++)
        {
            kineticEnergy(i,j) = 0.5*planets[j].M*pow(norm(vel.slice(i).col(j)),2);

            //potential energy from sun
            potentialEnergy(i,j) = -scale*planets[j].M/
            norm(pos.slice(i).col(j));
            //potential energy inbetween planets
            for(int k=j+1; k<numPlanets; k++)
            {
                double temp = -scale*planets[j].M*planets[k].M/
                norm(pos.slice(i).col(j) - pos.slice(i).col(k));
                potentialEnergy(i,j) += temp;
                potentialEnergy(i,k) += temp;
            }

            angularMomentum(i,j) = planets[j].M*
            norm(cross(pos.slice(i).col(j), vel.slice(i).col(j)));
        }
    }
}

void Verlet::coordinatesToFile(string filename)
{
    if (solved == false)
    {
        throw invalid_argument("Must run .solve() before .solveEnergy()");
    }
    ofstream myfile;
    myfile.open(filename);
    for(int j=0; j<N; j++)
    {
        for(int i=0; i<numPlanets; i++)
        {
            myfile << pos.slice(j).col(i)(0) << " "
                   << pos.slice(j).col(i)(1) << " "
                   << pos.slice(j).col(i)(2) << " "

                   << vel.slice(j).col(i)(0) << " "
                   << vel.slice(i).col(i)(1) << " "
                   << vel.slice(i).col(i)(2) << " ";
        }
        myfile << endl;
    }
    myfile.close();
}

void Verlet::energyToFile(string filename)
{
    if (solvedEnergy == false)
    {
        throw invalid_argument("Must run .solveEnergy() before .energyToFile()");
    }

    ofstream myfile;
    myfile.open(filename);
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<numPlanets; j++)
        {
            myfile << kineticEnergy(i,j) << " "
                   << potentialEnergy(i,j) << " "
                   << angularMomentum(i,j) << " ";
        }
        myfile << endl;
    }
    myfile.close();
}
