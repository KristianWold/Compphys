#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <vector>
using namespace std;
using namespace arma;

double scale = 39.485623;

vec newton(vec pos)
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

    Planet(){}
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
    vec t;
    cube pos;
    cube vel;
    vector<Planet> planets;
    int numPlanets;
    int N;

    void totalAcceleration(mat &totalAcc, int i)
    {
        totalAcc = zeros(3,numPlanets);
        for(int j=0; j<numPlanets; j++)
        {
            totalAcc.col(j) += newton(pos.slice(i).col(j));
            for(int k=j+1; k<numPlanets; k++)
            {
                mat temp = newton(pos.slice(i).col(j) - pos.slice(i).col(k));
                totalAcc.col(j) += planets[k].M*temp;
                totalAcc.col(k) -= planets[j].M*temp;
            }
        }
    }

public:
    Verlet(vector<Planet> p, int n)
    {
        planets = p;
        numPlanets = n;
    }

    void solve(vec acc(vec), double T, double dt)
    {
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
        totalAcceleration(totalAcc, 0);
        for(int i=0; i<N-1; i++)
        {
            pos.slice(i+1) = pos.slice(i) + vel.slice(i)*dt + 0.5*totalAcc*dt*dt;
            prevAcc = totalAcc;
            totalAcceleration(totalAcc, i+1);
            vel.slice(i+1) = vel.slice(i) + 0.5*(totalAcc + prevAcc)*dt;
        }
    }

    void writeToFile(string filename)
    {

        ofstream myfile;
        myfile.open(filename);
        for(int i=0; i<numPlanets; i++)
        {
            cout << i << endl;
            for(int j=0; j<N; j++)
            {
                myfile << pos.slice(j).col(i)(0) << " "
                       << pos.slice(j).col(i)(1) << " "
                       << pos.slice(j).col(i)(2) << " "

                       << vel.slice(j).col(i)(0) << " "
                       << vel.slice(i).col(i)(1) << " "
                       << vel.slice(i).col(i)(2) << endl;
            }
        }
    }
};

int main()
{
    vec pos = {1,0,0};
    vec vel = {0,2*M_PI,0};
    Planet Earth(pos, vel, 0.0001);
    vector<Planet> solarsystem = vector<Planet>{Earth};

    Verlet solver(solarsystem, 1);
    solver.solve(newton, 1, 0.001);
    solver.writeToFile("data.txt");
    return 0;
}
