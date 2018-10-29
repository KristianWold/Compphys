#include<stdlib.h>
#include<armadillo>
#include<fstream>
#include<cmath>
#include <map>
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
    int periodic(int x)
    {
        if(x == -1)
        {
            return L-1;
        }
        else if(x == L)
        {
            return 0;
        }
        else
        {
            return x;
        }
    }

public:
    int L;               //Dimenstionality LxL
    double J;            //Coupling constant
    double T;            //Temperature
    double energy;       //energy of the system
    double magnetization;//magnetization of the system
    double deltaE;       //change in energy when flipping single spin
    double deltaM;       //change in magnetization when flipping single spin
    Spins(){};
    Spins(int L, double T, double J)
    {
        this->L = L;
        this->T = T;
        this->J = J;
        acceptAmp.insert( pair<double,double>(-8,exp(-1/T*(-8))) );
        acceptAmp.insert( pair<double,double>(-4,exp(-1/T*(-4))) );
        acceptAmp.insert( pair<double,double>(0,1));
        acceptAmp.insert( pair<double,double>(4,exp(-1/T*(4))) );
        acceptAmp.insert( pair<double,double>(8,exp(-1/T*(8))) );
        ensemble = 2*randi<Mat<int>>(L,L,distr_param(0,1)) - ones<Mat<int>>(L,L);
        calcEnergy();
        calcMagnetization();
    }

    void calcEnergy()
    {
        energy = 0;
        for(int i=0; i<L-1; i++)
        {
            energy += -ensemble(L-1,i)*(ensemble(L-1,i+1) + ensemble(0,i));
            energy += -ensemble(i,L-1)*(ensemble(i+1,L-1) + ensemble(i,0));
            for(int j=0; j<L-1; j++)
            {
                energy += -ensemble(i,j)*(ensemble(i+1,j) + ensemble(i,j+1));
            }
        }
        energy += -ensemble(L-1,L-1)*(ensemble(L-1,0) + ensemble(0,L-1));
    }

    void calcMagnetization()
    {
        magnetization = 0;
        for(int i=0; i<L; i++)
        {
            for(int j=0; j<L; j++)
            {
                magnetization += ensemble(i,j);
            }
        }
    }

    void print()
    {
        cout << ensemble << endl;
    }

    void tryflip(double &aA)
    {
        x = randi<int>(distr_param(0,L-1));
        y = randi<int>(distr_param(0,L-1));
        deltaE = 2*ensemble(x,y)*(
                   ensemble(x,periodic(y+1)) + ensemble(x,periodic(y-1)) +
                   ensemble(periodic(x+1),y) + ensemble(periodic(x-1),y));
        deltaM = -2*ensemble(x,y);
        aA = acceptAmp.find(deltaE)->second; //Returns the amplitude associated
                                             //with the change in energy
    }
    void flip()
    {
        ensemble(x,y) *= -1;
    }
};

class MonteCarlo
{
private:
    Spins spins;
    double acceptAmp;

public:
    double energy;
    double av_energy;
    double e2;
    double magnetization;
    double av_magnetization;
    double Cv;

    MonteCarlo(Spins spins)
    {
        this->spins = spins;
    }
    void solve(int cycles)
    {
        ofstream myfile;
        myfile.open("data.txt");

        energy = spins.energy;
        av_energy = energy;
        e2 = energy*energy;

        magnetization = spins.magnetization;
        av_magnetization = magnetization;

        myfile << energy << endl;

        for(int i=1; i<cycles; i++)
        {
            //Sweeps over LxL spin matrix
            for(int j=0; j<spins.L*spins.L; j++)
            {
                spins.tryflip(acceptAmp);
                if(randu<double>() < acceptAmp)
                {
                    spins.flip();
                    energy += spins.deltaE;
                    magnetization += spins.deltaM;
                }
            }
            myfile << energy << "\n";

            av_energy += energy;
            e2 += energy*energy;

            av_magnetization += magnetization;
        }
        myfile.close();
        av_energy /= cycles;
        av_magnetization = abs(av_magnetization)/(cycles*spins.L*spins.L);
        e2 /= cycles;
        Cv = e2 - av_energy*av_energy;
    }
};

int main(int argc, char const *argv[])
{

    int cycles = atoi(argv[1]);
    int L = atoi(argv[2]);
    double T = atof(argv[3]);
    //double T;

    ofstream myfile;
    myfile.open("magnetization.txt");
    arma_rng::set_seed_random();

    Spins crystal(L,T,1);
    MonteCarlo MC(crystal);
    MC.solve(cycles);
    cout << MC.av_energy << endl;
    cout << MC.av_magnetization << endl;
    cout << MC.Cv << endl;

/*    for(int i = 0; i<=20; i++)
    {
        T = 2 + 0.3/20*i;
        Spins crystal(L,T,1);
        MonteCarlo MC(crystal);

        MC.solve(cycles);
        myfile << T << " " << MC.av_magnetization << "\n";
        cout << i << "/20" << endl;
    }
    myfile.close();
*/
    //system("python plot.py");

    return 0;
}
