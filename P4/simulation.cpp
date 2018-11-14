#include<stdlib.h>
#include<armadillo>
#include<fstream>
#include<cmath>
#include<map>
#include<random>
#include"classes/spins.h"
#include"classes/montecarlo.h"
#include"/usr/include/mpi/mpi.h"
using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
    int totalCycles = atoi(argv[1]);
    int cutoff = atoi(argv[2]);
    int L = atoi(argv[3]);
    double T = atof(argv[4]);
    int seed = atoi(argv[5]);
    int cycles;
    int* local;

    int numprocs, my_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


    cycles = totalCycles/numprocs + cutoff;

    mt19937_64 engine(seed + 1000*my_rank);
    //Spins crystal(L, T, 1, engine);                         //Disodered lattice
    Spins crystal(ones<Mat<int>>(L,L), L, T, 1, engine);  //Ordered lattice
    MonteCarlo MC(crystal);
    double start = MPI_Wtime();
    MC.solve(cycles, engine);
    double finish = MPI_Wtime();
    local = MC.energyAndMag;

    if(my_rank == 0)
    {
        MPI_Status status;
        ofstream file("results/data.dat", ofstream::binary);
        file.write(reinterpret_cast<const char*>(local),
                   2*cycles*sizeof(int));
        cout << "Numbers of accepted states: " << MC.accepted << endl;
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
        cout << finish - start << " seconds" << endl;
        cout << "Done!" << endl;
        ofstream meta;
        meta.open("results/meta.txt");
        meta << cycles << endl;
        meta << cutoff << endl;
        meta << numprocs << endl;
        meta << L << endl;
        meta << T << endl;
        meta << MC.accepted << endl;
        meta.close();
        system("python3 analyze.py");
    }
    return 0;
}
