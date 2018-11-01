#include <iostream>
#include <map>
#include <armadillo>
#include<fstream>
#include <random>
#include"/usr/include/mpi/mpi.h"
using namespace std;
using namespace arma;

int *intvector(int row)
{
  int * A;
  A = new int[row];
  return (int *)A;
}

int main(int nargs, char* args[])
{
    int* local;

    int numprocs, my_rank;
    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    for(int i=1; i<10; i++)
    {
        
    }

    if(my_rank == 0)
    {
        MPI_Status status;
        cout << local << endl;
        for(int j=1; j<numprocs; j++)
        {

            MPI_Recv(&local, 1, MPI_INT, MPI_ANY_SOURCE, 500,
                     MPI_COMM_WORLD, &status);
            cout << local << endl;
        }
    }
    else
    {
        MPI_Send(&local, 1, MPI_INT, 0, 500, MPI_COMM_WORLD);
    }
    MPI_Finalize();

    return 0;
}
