#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <math.h>
#include <string>
using namespace std;
using namespace arma;

double scale = 4*M_PI*M_PI/(365.25*365.25);

vec acceleration(vec pos)
{
    double M = 1;
    double rCube = pow(norm(pos), 3);
    return -scale*M/rCube*pos;
}

void euler(vec pos, vec vel, vec acc(vec), int N, double dt, int sampleN)
{
    ofstream myfile;
    myfile.open("twoBody.txt");
    for(int i=0; i<N-1; i++)
    {
        vel = vel + acc(pos)*dt;
        pos = pos + vel*dt;
        if(i%sampleN == 0)
        {
            myfile << pos(0) << " " << pos(1) << " " << pos(2);
            myfile << vel(0) << " " << vel(1) << " " << vel(2);
            myfile << "\n";
        }
    }
    myfile.close();
}

void velocityVerlet(vec pos, vec vel, vec acc(vec), int N, double dt, int sampleN)
{
    ofstream myfile;
    myfile.open("twoBody.txt");

    vec prevPos(3,fill::zeros);
    for(int i=0; i<N-1; i++)
    {
        prevPos = pos;
        pos = pos + vel*dt + 0.5*acc(pos)*dt*dt;
        vel = vel + 0.5*(acc(pos) + acc(prevPos))*dt;
        if(i%sampleN == 0)
        {
            myfile << pos(0) << " " << pos(1) << " " << pos(2);

            myfile << "\n";
        }
    }
    myfile.close();
}

int main(int argc, char const *argv[])
{
    int solveType = atoi(argv[1]);  //1 for forward Euler, 2 for velocityVerlet
    double T = atof(argv[2]);       //final time
    int N = atoi(argv[3]);          //Number of integration points
    int sampleN = atoi(argv[4]);    //Increments to write to file

    double dt = double(T/N);        //timestep

    vec pos({1,0,0});
    vec vel({0,2*M_PI/365.25,0});
    if (solveType == 1)
    {
        euler(pos, vel, acceleration, N, dt, sampleN);
    }
    if (solveType == 2)
    {
        velocityVerlet(pos, vel, acceleration, N, dt, sampleN);
    }
}
