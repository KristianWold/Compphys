#include <iostream>
#include <chrono>
#include <math.h>
#include <armadillo>
#include <fstream>

#include "func.h"

using namespace arma;
using namespace std;
using namespace chrono;

vec harmOsc(vec x, double w)
{
	return w*w*x%x;
}

vec coulomb(vec x, double w)
{
	return w*w*x%x + 1/x;
}

int main(int argc, char *argv[])
{
	if (argc < 5)
	{
		cout << "Bad usage" << endl;
		cout << "Please supply:" << endl;
		cout << "	number of discrete point n" << endl;
		cout << "	starting point p0" << endl;
		cout << "	end point pN" << endl;
		cout << "	oscillator frequencys w" << endl;
		return 1;
	}
	int n = atof(argv[1]);		//Grid-points
	double p0 = atof(argv[2]);	//Start-point
	double pN = atof(argv[3]);	//End-points
	double w = atof(argv[4]);	//Oscillator frequency

	double h = (pN - p0)/(n+1);			//Step-size
	vec x = linspace(p0+h, pN-h, n);	//Spacing

	vec eigval;
	mat A;
	mat eigvec1, eigvec2;

	//non-interacting case
	A = makeMatrix(harmOsc(x, w), h, n);
	solveJacobi(A, eigval, eigvec1, n);

	//interacting case
	A = makeMatrix(coulomb(x, w), h, n);
	solveJacobi(A, eigval, eigvec2, n);

	ofstream myfile;
	myfile.open("eigenvec.txt");
	for(int i=0; i<n; i++)
	{
		myfile << x(i) <<
		" " << eigvec1(i,0)/sqrt(h)<<
		" " << eigvec2(i,0)/sqrt(h)<< endl;
	}
	myfile.close();

	string command = string("python plot.py eigenvectors ") + string(argv[4]);
	system(command.c_str());

	return 0;
}
