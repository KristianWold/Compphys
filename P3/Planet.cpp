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
        for(int j=0; j<numPlanets-1; j++)
        {
            totalAcc[j] += newton(pos.slice(j).col(i));
            for(int k=j+1; k<numPlanets; k++)
            {
                mat temp = newton(pos.slice(j).col(i) - pos.slice(k).col(i));
                totalAcc[j] += planets[k].M*temp;
                totalAcc[k] -= planets[j].M*temp;
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
        pos = zeros(3, n, N);
        vel = zeros(3, n, N);
        for(int i=0; i<numPlanets; i++)
        {
            pos.slice(0).col(i) = planets[i].pos;
            vel.slice(0).col(i) = planets[i].vel;
        }

        mat totalAcc(N, n, fill::zeros);
        mat prevAcc(N, n, fill::zeros);
        totalAcceleration(totalAcc, 0);
        for(int i=0; i<N-1; i++)
        {
            pos.col(i+1) = pos.col(i) + vel.col(i)*dt + 0.5*acc(pos.col(i))*dt*dt;
            vel.col(i+1) = vel.col(i) + 0.5*(acc(pos.col(i+1)) + acc(pos.col(i)))*dt;
        }


    }
/*
    void writeToFile(string filename)
    {

        ofstream myfile;
        myfile.open(filename);
        for(int i=0; i<N; i++)
        {
            myfile << pos.col(i)(0) << " " << pos.col(i)(1) << " " << pos.col(i)(2);
            myfile << vel.col(i)(0) << " " << vel.col(i)(1) << " " << vel.col(i)(2);
            myfile << endl;
        }
    }*/
};

int main()
{
    vec pos = {1,0,0};
    vec vel = {0,2*M_PI,0};
    Planet Earth(pos, vel);
    vector<Planet> solarsystem(Earth);

    Verlet solver(solarsystem);
    solver.solve(newton, 1, 0.001);
    solver.writeToFile("data.txt");
    return 0;
}
