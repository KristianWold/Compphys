#include <iostream>
#include <math.h>
#include <fstream>
#include "func.cpp"


int main(int argc, char *argv[])
{

    //meta.bin tells python the dimmensions of the vectors
    int n[argc];
    n[0] = argc-1;
    for(int i = 1; i<argc; i++)
    {
        n[i] = atoi(argv[i]);
    }
    FILE * pFile;
    pFile = fopen ("meta.bin", "wb");
    fwrite (n, sizeof(int), sizeof(n), pFile);
    fclose (pFile);

    //writes data to binary file so python can read it
    pFile = fopen ("myfile.bin", "wb");
    for (int i = 1; i<=n[0]; i++)
    {
        fwrite (solveB(n[i], source), sizeof(double), n[i], pFile);
        //fwrite (solveC(n[i], source), sizeof(double), n[i], pFile);
    }
    fclose (pFile);
    return 0;
}
