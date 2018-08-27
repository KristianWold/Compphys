#include <iostream>
#include <math.h>
#include <fstream>
#include "func.cpp"


int main()
{

  //meta.bin tells python the dimmensions of the vectors
  int n[] = {3,12,102,1002};
  FILE * pFile;
  pFile = fopen ("meta.bin", "wb");
  fwrite (n, sizeof(int), sizeof(n), pFile);
  fclose (pFile);

  //writes data to binary file so python can read it
  pFile = fopen ("myfile.bin", "wb");
  for (int i = 1; i<n[0]+1; i++)
  {
    fwrite (solveC(n[i]), sizeof(double), n[i], pFile);
  }
  fclose (pFile);

  solveB(1e3);
  solveC(1e3);
  return 0;
}
