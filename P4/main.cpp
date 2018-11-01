#include<stdlib.h>
#include<armadillo>
#include<fstream>
#include<cmath>
#include<map>
#include<random>
#include"/usr/include/mpi/mpi.h"
using namespace std;
using namespace arma;

int *intvector(int row)
{
  int * A;
  A = new int[row];
  return (int *)A;
}

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
    void init(int L, double T, double J, mt19937 &engine)
    {
        this->L = L;
        this->T = T;
        this->J = J;

        rand_spin  = uniform_int_distribution<int>(0,1);
        rand_coord = uniform_int_distribution<int>(0,L-1);

        acceptAmp.insert( pair<double,double>(-8,exp(-1/T*(-8))) );
        acceptAmp.insert( pair<double,double>(-4,exp(-1/T*(-4))) );
        acceptAmp.insert( pair<double,double>(0,1));
        acceptAmp.insert( pair<double,double>(4,exp(-1/T*(4))) );
        acceptAmp.insert( pair<double,double>(8,exp(-1/T*(8))) );
    }

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


    Spins(){};
    Spins(int L, double T, double J, mt19937 &engine)
    {
        init(L, T, J, engine);
        ensemble = zeros<Mat<int>>(L,L);
        for(int i=0; i<L; i++)
        {
            for(int j=0; j<L; j++)
            {
                ensemble(i,j) = 2*rand_spin(engine) - 1;
            }
        }
        calcEnergy();
        calcMagnetization();
    }

    Spins(Mat<int> ensemble, int L, double T, double J, mt19937 &engine)
    {
        init(L, T, J, engine);
        this->ensemble = ensemble;
        calcEnergy();
        calcMagnetization();
    }

    void calcEnergy()
    {
        energy = 0;
        for(int i=0; i<L-1; i++)
        {
            energy -= ensemble(L-1,i)*(ensemble(L-1,i+1) + ensemble(0,i));
            energy -= ensemble(i,L-1)*(ensemble(i+1,L-1) + ensemble(i,0));
            for(int j=0; j<L-1; j++)
            {
                energy -= ensemble(i,j)*(ensemble(i+1,j) + ensemble(i,j+1));
            }
        }
        energy -= ensemble(L-1,L-1)*(ensemble(L-1,0) + ensemble(0,L-1));
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

    void tryflip(double &aA, mt19937 &engine)
    {
        x = rand_coord(engine);
        y = rand_coord(engine);
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
    int* energyAndMag;

    uniform_real_distribution<float> rand_float;

    MonteCarlo(Spins spins)
    {
        this->spins = spins;
        rand_float = uniform_real_distribution<float>(0,1);

    }
    void solve(int cycles, mt19937 &engine)
    {
        energyAndMag = intvector(2*cycles);

        energyAndMag[0]=spins.energy;
        energyAndMag[cycles] = spins.magnetization;

        for(int i=1; i<cycles; i++)
        {
            //Sweeps over LxL spin matrix
            energyAndMag[i] = energyAndMag[i-1];
            energyAndMag[i+cycles] = energyAndMag[i+cycles-1];
            for(int j=0; j<spins.L*spins.L; j++)
            {
                spins.tryflip(acceptAmp, engine);
                if(rand_float(engine) < acceptAmp)
                {
                    spins.flip();
                    energyAndMag[i] += spins.deltaE;
                    energyAndMag[i+cycles] += spins.deltaM;
                }
            }
        }
        cout << "Done!" << endl;
    }
};

int main(int argc, char *argv[])
{
    int cycles = atoi(argv[1]);
    int L = atoi(argv[2]);
    double T = atof(argv[3]);
    int* local;

    int numprocs, my_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    mt19937 engine(my_rank+10);
    Spins crystal(L, T, 1, engine);
    MonteCarlo MC(crystal);
    MC.solve(cycles, engine);
    local = MC.energyAndMag;

    if(my_rank == 0)
    {
        MPI_Status status;
        ofstream file("results/data.dat", ofstream::binary);
        file.write(reinterpret_cast<const char*>(local),
                   2*cycles*sizeof(int));

        for(int i=1; i<numprocs; i++)
        {
            MPI_Recv(local, 2*cycles, MPI_INT, MPI_ANY_SOURCE, 500,
                     MPI_COMM_WORLD, &status);
            file.write(reinterpret_cast<const char*>(local),
                       2*cycles*sizeof(int));
        }
        file.close();
    }
    else
    {
        MPI_Send(local, 2*cycles, MPI_INT, 0, 500, MPI_COMM_WORLD);
    }
    MPI_Finalize();

    if(my_rank == 0)
    {
        ofstream meta;
        meta.open("results/meta.txt");
        meta << cycles << endl;
        meta << numprocs << endl;
        meta << L << endl;
        meta << T << endl;
        meta.close();
        system("python3 analyze.py");
    }
    return 0;
}
