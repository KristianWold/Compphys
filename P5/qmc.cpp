#include<cmath>
#include<random>
#include<iostream>
#include<iomanip>

using namespace std;

double acceptAmp(double alpha, double omega, double* delta, double* r)
{
    double expo = delta[0]*(delta[0]+2*r[0]) + delta[1]*(delta[1]+2*r[1]) + delta[2]*(delta[2]+2*r[2]);
    return exp(-alpha*omega*expo/2);
}

int main(int argc, char const *argv[])
{
    int N = 2;
    int dim = 3;
    double alpha = 2;
    double omega = 1;
    double step = 3/sqrt(alpha*omega);
    int accepted = 0;
    double E = 0;
    double R = 0;

    double* delta = new double[dim];
    double** r = new double*[2];
    r[0] = new double[3];
    r[1] = new double[3];

    for (int i=0; i<N; i++)
    {
        for (int j=0; j<dim; j++)
        {
            r[i][j] = 0;
        }
    }

    mt19937_64 engine;
    uniform_real_distribution<double> myRandu(0,1);

    int M = 100000;
    for(int i=0; i<M; i++)
    {
        for(int j=0; j<N; j++)
        {
            delta[0] = (myRandu(engine) - 0.5)*step;
            delta[1] = (myRandu(engine) - 0.5)*step;
            delta[2] = (myRandu(engine) - 0.5)*step;

            if(myRandu(engine) < acceptAmp(alpha, omega, delta, r[j]))
            {
                accepted++;
                r[j][0] += delta[0];
                r[j][1] += delta[1];
                r[j][2] += delta[2];
            }
        }
        R = 0;
        for (int i=0; i<N; i++)
        {
            for (int j=0; j<dim; j++)
            {
                R += r[i][j]*r[i][j];
            }
        }

        double R12 = sqrt(abs(r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[0][2]*r[0][2]
                            - r[1][0]*r[1][0] + r[1][1]*r[1][1] + r[1][2]*r[0][2]));

        E+= 0.5*omega*omega*R*(1-alpha*alpha) + 3*alpha*omega;
    }
    E /= M;

    cout << accepted << endl;
    cout <<setprecision(10) << E << endl;

    return 0;
}
