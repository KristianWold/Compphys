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

Verlet::Verlet(vector<Planet> p, int n, double s)
{
    planets = p;
    numPlanets = n;
    scale = s;
}

void Verlet::solve(vec acc(vec, vec), double T, int N, int sampleN)
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
        pos.col(i) = planets[i].pos;
        vel.col(i) = planets[i].vel;
    }

    ofstream myfile;
    myfile.open("data.txt");

    mat totalAcc(3, numPlanets, fill::zeros);
    mat prevAcc(3, numPlanets, fill::zeros);
    totalAcceleration(totalAcc, acc, 0);
    for(int i=0; i<N-1; i++)
    {
        pos = pos + vel*dt + 0.5*totalAcc*dt*dt;
        prevAcc = totalAcc;
        totalAcceleration(totalAcc, acc, i+1);
        vel = vel + 0.5*(totalAcc + prevAcc)*dt;

        if (i%sampleN == 0)
        {
            for(int j=0; j<numPlanets; j++)
            {
                myfile << pos.col(j)(0) << " "
                       << pos.col(j)(1) << " "
                       << pos.col(j)(2) << " "

                       << vel.col(j)(0) << " "
                       << vel.col(j)(1) << " "
                       << vel.col(j)(2) << " ";

                kineticEnergy(countSample,j) = 0.5*planets[j].M*pow(norm(vel.col(j)),2);

                //potential energy from sun
                potentialEnergy(countSample,j) = -scale*planets[j].M/
                norm(pos.col(j));
                //potential energy inbetween planets
                for(int k=j+1; k<numPlanets; k++)
                {
                    double temp = -scale*planets[j].M*planets[k].M/
                    norm(pos.col(j) - pos.col(k));
                    potentialEnergy(countSample,j) += temp;
                    potentialEnergy(countSample,k) += temp;
                }

                angularMomentum(countSample,j) = planets[j].M*
                norm(cross(pos.col(j), vel.col(j)));
            }
            myfile << "\n";

            countSample++;

        }
    }
    myfile.close();
}
/*
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
*/
