#include <iostream>
#include <iomanip>
#include <chrono>
#include <math.h>
#include <armadillo>
#include <fstream>

#include "jacobi.h"

using namespace arma;
using namespace std;
using namespace chrono;

//=============================================================================
vec harmOsc(vec x, double w)
//Describes one electron(or two non-interaction electons)
//in a harmonic oscillator
//----------------------------------------------------------------------------
//x - interval the potential is defined
//w - oscillator frequency
{
	return w*w*x%x;
}

vec coulomb(vec x, double w)
//Describes twi interaction electrons in a harmonic
//----------------------------------------------------------------------------
//x - interval the potential is defined
//w - oscillator frequency
{
	return w*w*x%x + 1/x;
}

int main(int argc, char const *argv[])
{
	//Checks if the user is supplying enough commandline arguments
	if (argc < 6)
	{
		cout << "Bad usage" << endl;
		cout << "Please supply:" << endl;
		cout << "	number of discrete point n" << endl;
		cout << "	starting point p0" << endl;
		cout << "	end point pN" << endl;
		cout << "	oscillator frequency w" << endl;
		cout << "	some potential: 1 for harmOsc, 2 for coulomb" << endl;
		return 1;
	}
	int n = atof(argv[1]);		//Grid-points
	double p0 = atof(argv[2]);	//Start-point
	double pN = atof(argv[3]);	//End-points
	double w = atof(argv[4]);	//Oscillator frequency
	int typePotential = atoi(argv[5]);

	double h = (pN - p0)/(n+1);			//Step-size
	vec x = linspace(p0+h, pN-h, n);	//Spacing

	vec potential;
	if (typePotential == 1)
	{
		potential = harmOsc(x,w);
	}
	else if (typePotential == 2)
	{
		potential = coulomb(x,w);
	}
	else
	{
		cout << "Bad usage" << endl;
		cout << "Please use 1 for harmOsc, 2 for coulomb" << endl;
		return 1;
	}

	mat A;
	vec eigval;
	mat eigvec;

	A = makeMatrix(potential, h, n);
	//Calls the method to yield eigenvalues and eigenvectors
	solveJacobi(A, eigval, eigvec, n);

	cout << setprecision(7) << "E0: " << eigval(0) << endl;
	cout << setprecision(7) << "E1: " << eigval(1) << endl;
	cout << setprecision(7) << "E2: " << eigval(2) << endl;
	cout << setprecision(7) << "E3: " << eigval(3) << endl;
	return 0;
}
