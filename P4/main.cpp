#include<stdlib.h>
#include<armadillo>
#include<fstream>
#include<cmath>
using namespace std;
using namespace arma;


double k = 1.38064852*1e-23;

class Spins
{
private:
    Mat<int> ensemble;
    int L;
    double J;

public:
    Spins(){};
    Spins(int L, double J)
    {
        this->L = L;
        this->J = J;
        newEnsemble();
    }
    double energy()
    {
        double E = 0;
        for(int i=0; i<L-1; i++)
        {
            E += ensemble(L-1,i)*(ensemble(L-1,i+1) + ensemble(0,i));
            E += ensemble(i,L-1)*(ensemble(i+1,L-1) + ensemble(i,0));
            for(int j=0; j<L-1; j++)
            {
                E += ensemble(i,j)*(ensemble(i+1,j) + ensemble(i,j+1));
            }
        }
        E += ensemble(L-1,L-1)*(ensemble(L-1,0) + ensemble(0,L-1));

        return J*E;
    }

    void print()
    {
        cout << ensemble << endl;
    }
    void newEnsemble()
    {
        ensemble = 2*randi<Mat<int>>(L, L, distr_param(0,1)) - ones<Mat<int>>(L, L);
    }
};

class MonteCarlo
{
private:
    Spins spins;
    bool accept(double energy, double newEnergy, double T)
    {
        if(newEnergy <= energy)
        {
            return true;
        }

        else
        {
            if (randu<double>() < exp(1/T*(energy - newEnergy)))
            {
                return true;
            }
        }
        return false;
    }
public:
    double energy;
    double av_energy;
    double newEnergy;
    MonteCarlo(Spins spins)
    {
        this->spins = spins;
    }
    void solve(int N, double T)
    {
        ofstream myfile;
        myfile.open("data.txt");
        energy = spins.energy();
        av_energy = energy;
        for(int i=1; i<N; i++)
        {
            spins.newEnsemble();
            newEnergy = spins.energy();
            while(!accept(energy, newEnergy, T))
            {
                spins.newEnsemble();
                newEnergy = spins.energy();
            }
            energy = newEnergy;
            av_energy += energy;
            myfile << energy << "\n";
        }
        myfile.close();
        av_energy /= N;
    }
};

int main(int argc, char const *argv[])
{
    int N = atoi(argv[1]);
    int L = atoi(argv[2]);
    double T = atof(argv[3]);
    arma_rng::set_seed_random();
    Spins crystal(L,1);
    MonteCarlo MC(crystal);
    MC.solve(N,T);
    cout << MC.av_energy << endl;

    return 0;
}
