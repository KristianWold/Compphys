#include <iostream>
#include <iomanip>
#include <armadillo>
using namespace std;
using namespace arma;

double G = 1.;
double M = 1.;

vec acceleration(vec pos)
{
    double rCube = pow(norm(pos), 3);
    return G*M/rCube*pos;
}

void Euler(mat &pos, mat &vel, vec acc(vec), int N, double dt)
{
    for(int i=0; i<N; i++)
    {
        vel.col(i+1) = vel.col(i) + acc(pos.col(i))*dt;
        pos.col(i+1) = pos.col(i) + vel.col(i)*dt;
    }
}

int main()
{
    double T = 10;
    double dt = 0.01;
    int N = int(T/dt);

    mat pos(3,N,fill::zeros); pos.col(0) =
    {1,2,3};
    //mat vel(3,N,fill::zeros); vel.col(0) = {1,2,3};

    //Euler(pos, vel, acceleration, T, dt)

    return 0;
}
