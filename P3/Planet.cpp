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
    vec t;
    cube pos;
    cube vel;
    vector<Planet> planets;
    int numPlanets;
    int N;

    void totalAcceleration(mat &totalAcc, int i)
    {
        totalAcc = zeros(3, numPlanets);
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
};

int main(int argc, char const *argv[])
{
    Planet Earth(vec({1,0,0}), vec({0,2*M_PI/365.25 ,0}),1./333000);
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
    solver.writeToFile("data.txt");

    system("python plot.py");

    return 0;
}
