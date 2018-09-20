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
		int n = 500;
		double pN = 15;
		double p0 = 0.0;
		double w = 0.1;

		double h = (pN - p0)/(n+1);
		double h2 = h*h;
		vec x = linspace(p0+h, pN-h, n);

		//Col<double> V = w*w*x%x;// +1/x;

		auto start = high_resolution_clock::now();
		mat eigvec1 = SolveArma(n, h, x, HarmonicOsc, w);
		mat eigvec2 = SolveArma(n, h, x, Interacting, w);
		auto finish = high_resolution_clock::now();

		cout << duration<double>(finish - start).count() << endl;

		ofstream myfile;
		myfile.open("eigenvec.txt");
		for(int i=0; i<n; i++)
		{
				myfile << x(i) << " " << eigvec1(i,0)/sqrt(h) << " " << eigvec2(i,0)/sqrt(h)<< endl;
		}
		myfile.close();

		return 0;
}
