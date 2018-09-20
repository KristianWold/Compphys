#include <iostream>
#include <chrono>
#include <math.h>
#include <armadillo>
#include <fstream>
#include <chrono>
#include "func.cpp"
using namespace arma;
using namespace std;
using namespace chrono;

inline vec HarmonicOsc(double w,vec x)
{
		return w*w*x%x;
}

inline vec Interacting(double w,vec x)
{
		return w*w*x%x + 1/x;
}


int main()
{
		int n = 200;
		double pN = 3;
		double p0 = 0.0;

		double h = (pN - p0)/(n+1);
		vec x = linspace(p0+h, pN-h, n);

		vec eig = SolveJacobi(h, x, HarmonicOsc, 1, n);

		cout << eig(0) << endl;
		cout << eig(1) << endl;
		cout << eig(2) << endl;

/*
		n = 700;
		pN = 20;
	 	p0 = 0.0;

		h = (pN - p0)/(n+1);
		h2 = h*h;
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
