#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
using namespace std;
using namespace arma;

double scale = 39.485623;

vec acceleration(vec pos)
{
    double M = 1;
    double rCube = pow(norm(pos), 3);
    return -scale*M/rCube*pos;
}

void euler(mat &pos, mat &vel, vec acc(vec), int N, double dt)
{
    for(int i=0; i<N-1; i++)
    {
        vel.col(i+1) = vel.col(i) + acc(pos.col(i))*dt;
        pos.col(i+1) = pos.col(i) + vel.col(i)*dt;
    }
}

void velocityVerlet(mat &pos, mat &vel, vec acc(vec), int N, double dt)
{
    for(int i=0; i<N-1; i++)
    {
        pos.col(i+1) = pos.col(i) + vel.col(i)*dt + 0.5*acc(pos.col(i))*dt*dt;
        vel.col(i+1) = vel.col(i) + 0.5*(acc(pos.col(i+1)) + acc(pos.col(i)))*dt;
    }
}

int main(int argc, char const *argv[])
{
    int solveType = atoi(argv[1]);
    double T = atof(argv[2]);
    double dt = atof(argv[3]);

    int N = int(T/dt);

    mat pos(3,N,fill::zeros); pos.col(0) = vec({1,0,0});
    mat vel(3,N,fill::zeros); vel.col(0) = vec({0,2*M_PI,0});
    if (solveType == 1)
    {
        euler(pos, vel, acceleration, N, dt);
    }
    if (solveType == 2)
    {
        velocityVerlet(pos, vel, acceleration, N, dt);
    }


    ofstream myfile;
    myfile.open("data.txt");
    for(int i=0; i<N; i++)
    {
        myfile << pos.col(i)(0) << " " << pos.col(i)(1) << " " << pos.col(i)(2);
        myfile << vel.col(i)(0) << " " << vel.col(i)(1) << " " << vel.col(i)(2);
        myfile << endl;
    }
    system("python plot.py");
    return 0;
}
