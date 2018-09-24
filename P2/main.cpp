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

	int n = 200;
	double p0 = 0.0;
	double pN = 5;

	double h = (pN - p0)/(n+1);
	vec x = linspace(p0+h, pN-h, n);

	mat A;
	vec eigval;
	mat eigvec;
/*
	A = makeMatrix(x, 1, harmonicOsc, h, n);
	solveJacobi(A, eigval, eigvec, n);

	cout << eigval(0) << endl;
	cout << eigval(1) << endl;
	cout << eigval(2) << endl;
*/

	n = 100;
	p0 = 0.0;
	pN = 15;

	h = (pN - p0)/(n+1);
	x = linspace(p0+h, pN-h, n);

	mat eigvec1, eigvec2, eigvec3, eigvec4;

	A = makeMatrix(x, 1, harmonicOsc, h, n);
	solveJacobi(A, eigval, eigvec1, n);

	A = makeMatrix(x, 1, interacting, h, n);
	solveJacobi(A, eigval, eigvec2, n);

	A = makeMatrix(x, 0.3, interacting, h, n);
	solveJacobi(A, eigval, eigvec3, n);

	A = makeMatrix(x, 0.06, interacting, h, n);
	solveJacobi(A, eigval, eigvec4, n);

	ofstream myfile;
	myfile.open("eigenvec.txt");
	for(int i=0; i<n; i++)
	{
			myfile << x(i) << " " << eigvec1(i,0)/sqrt(h) << " " << eigvec2(i,0)/sqrt(h)<<
			        " " << eigvec3(i,0)/sqrt(h) << " " << eigvec4(i,0)/sqrt(h)<< endl;
	}
	myfile.close();

	return 0;
}
