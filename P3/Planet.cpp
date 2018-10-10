#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
using namespace std;
using namespace arma;

double scale = 4*M_PI*M_PI/(365.25*365.25);

vec newton(vec pos, vec vel)
{
    double rCube = pow(norm(pos), 3);
    return -scale/rCube*pos;
}

class Planet
{
public:
    vec pos;
    vec vel;
    double M;

    Planet()
    {
        pos = zeros(3);
        vel = zeros(3);
        M = 0;
    }

    Planet(vec position, vec velocity, double mass)
    {
        pos = position;
        vel = velocity;
        M = mass;
    }
};

class Verlet
{
private:
    //keeps track of if the methods solve() and solveEnergy() have been called
    bool solved = false;
    bool solvedEnergy = false;

    vec t;
    cube pos;
    cube vel;
    vector<Planet> planets;
    int numPlanets;
    int N;

    void totalAcceleration(mat &totalAcc, vec acc(vec, vec), int i)
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

public:
    mat kineticEnergy;
    mat potentialEnergy;
    mat angularMomentum;

    Verlet(vector<Planet> p, int n)
    {
        planets = p;
        numPlanets = n;
    }

    void solve(vec acc(vec, vec), double T, double dt)
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

    void solveEnergy()
    {
        if (solved == false)
        {
            throw invalid_argument("Must run .solve() first");
        }
        solvedEnergy = true;
        kineticEnergy = zeros(numPlanets, N);
        potentialEnergy = zeros(numPlanets, N);
        angularMomentum = zeros(numPlanets, N);

        for(int i=0; i<N; i++)
        {
            for(int j=0; j<numPlanets; j++)
            {
                kineticEnergy(j,i) = 0.5*planets[j].M*pow(norm(vel.slice(i).col(j)),2);

                //potential energy from sun
                potentialEnergy(j,i) = -scale*planets[j].M/
                norm(pos.slice(i).col(j));
                //potential energy inbetween planets
                for(int k=j+1; k<numPlanets; k++)
                {
                    double temp = -scale*planets[j].M*planets[k].M/
                    norm(vel.slice(i).col(j) - vel.slice(i).col(k));
                    potentialEnergy(j,i) += temp;
                    potentialEnergy(k,i) += temp;
                }

                angularMomentum(j,i) = planets[j].M*
                norm(cross(pos.slice(i).col(j), vel.slice(i).col(j)));
            }
        }
    }

    void coordinatesToFile(string filename)
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

    void energyToFile(string filename)
    {
        if (solvedEnergy == false)
        {
            throw invalid_argument( "Must run .solveEnergy() before .energyToFile()" );
        }

        ofstream myfile;
        myfile.open(filename);
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<numPlanets; j++)
            {
                myfile << kineticEnergy(j,i) << " "
                       << potentialEnergy(j,i) << " "
                       << angularMomentum(j,i) << " ";
            }
            myfile << endl;
        }
        myfile.close();
    }
};

int main(int argc, char const *argv[])
{
    Planet Earth(vec({1,0,0}), vec({0,M_PI/365.25,0}),1./333000);
    //Planet Jupiter; Jupiter.M = 0.0009543;

    vector<Planet> solarsystem = vector<Planet>{Earth};

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
    Verlet solver(solarsystem, 1);
    solver.solve(newton, atof(argv[1]), atof(argv[2]));
    solver.solveEnergy();

    solver.coordinatesToFile("coordinates.txt");
    solver.energyToFile("energy.txt");
    int a = system("python plot.py");

    return 0;
}
