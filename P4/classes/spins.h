#ifndef SPINS_H
#define SPINS_H

#include<stdlib.h>
#include<armadillo>
#include<fstream>
#include<cmath>
#include<map>
#include<random>
using namespace std;
using namespace arma;

class Spins
{
private:
    Mat<int> ensemble;  //Matrix of the spins
    int x, y;           //Coordinates of chosen spin to be flipped
    map<double, double> acceptAmp; //Map that relates change of energy to
                                   //the appropirate acceptance amplitude

    //takes care of the periodic boundary conditions
    int periodic(int x);

    void init(int L, double T, double J, mt19937_64 &engine);
public:
    int L;               //Dimenstionality LxL
    double J;            //Coupling constant
    double T;            //Temperature
    double energy;       //energy of the system
    double magnetization;//magnetization of the system
    double deltaE;       //change in energy when flipping single spin
    double deltaM;       //change in magnetization when flipping single spin

    uniform_int_distribution<int> rand_spin;
    uniform_int_distribution<int> rand_coord;


    Spins();
    Spins(int L, double T, double J, mt19937_64 &engine);
    Spins(Mat<int> ensemble, int L, double T, double J, mt19937_64 &engine);

    void calcEnergy();

    void calcMagnetization();

    void print();

    void tryflip(double &aA, mt19937_64 &engine);
    void flip();
};

#endif
