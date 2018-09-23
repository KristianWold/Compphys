#include <iostream>
#include <chrono>
#include <math.h>
#include <armadillo>
#include <fstream>

#include "func.h"

using namespace arma;
using namespace std;
using namespace chrono;

vec harmonicOsc(vec x, double w)
{
	return w*w*x%x;
}

vec interacting(vec x, double w)
{
	return w*w*x%x + 1/x;
}

int main()
{
	int n = 300;
	double pN = 5;
	double p0 = 0.0;

	double h = (pN - p0)/(n+1);
	vec x = linspace(p0+h, pN-h, n);

	mat A = makeMatrix(x, 1, harmonicOsc, h, n);

	vec eig = solveJacobi(A, n);

	cout << eig(0) << endl;
	cout << eig(1) << endl;
	cout << eig(2) << endl;

/*
		n = 700;
		pN = 20;
		p0 = 0.0;

		h = (pN - p0)/(n+1);
		x = linspace(p0+h, pN-h, n);

		mat eigvec1 = SolveArma(n, h, x, HarmonicOsc, 1);
		mat eigvec2 = SolveArma(n, h, x, Interacting, 1);
		mat eigvec3 = SolveArma(n, h, x, Interacting, 0.3);
		mat eigvec4 = SolveArma(n, h, x, Interacting, 0.05);

		ofstream myfile;
		myfile.open("eigenvec.txt");
		for(int i=0; i<n; i++)
		{
				myfile << x(i) << " " << eigvec1(i,0)/sqrt(h) << " " << eigvec2(i,0)/sqrt(h)<<
				        " " << eigvec3(i,0)/sqrt(h) << " " << eigvec4(i,0)/sqrt(h)<< endl;
		}
		myfile.close();
*/
		return 0;
}
