#include<stdlib.h>
#include<armadillo>
#include<fstream>
#include<cmath>
#include<map>
#include<random>
#include"spins.h"
#include"montecarlo.h"

using namespace std;
using namespace arma;

MonteCarlo::MonteCarlo(Spins spins)
{
    this->spins = spins;
    rand_float = uniform_real_distribution<float>(0,1);
}
void MonteCarlo::solve(int cycles, mt19937 &engine)
{
    energyAndMag = new int[2*cycles];

    energyAndMag[0] = E = spins.energy;
    energyAndMag[cycles] = M = abs(spins.magnetization);

    for(int i=1; i<cycles; i++)
    {
        //Sweeps over LxL spin matrix
        for(int j=0; j<spins.L*spins.L; j++)
        {
            spins.tryflip(acceptAmp, engine);
            if(rand_float(engine) < acceptAmp)
            {
                spins.flip();
                E += spins.deltaE;
                M += spins.deltaM;
            }
        }
        energyAndMag[i] = E;
        energyAndMag[i+cycles] = abs(M);
    }
}
