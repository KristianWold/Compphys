//#include <func.h>
#include <iostream>
#include <math.h>
#include <armadillo>
#include <utility>
#include "func.h"

using namespace arma;
using namespace std;

mat makeMatrix(vec x, double w, vec V(vec x, double w), double h, int n)
{
	double h2 = h*h;
	mat A(n,n, fill::zeros);
	A.diag(-1).fill(-1/h2);
	A.diag(0).fill(2/h2);
	A.diag(0) += V(x,w);
	A.diag(1).fill(-1/h2);

	return A;
}

//Finds the largest off-diagonal element above the diagonal
void max_element(mat &A, int &k, int &l, double &max, int n)
{
	max = 0;
	for(int i = 0; i<n; i++)
	{
		for(int j = i+1; j<n; j++)
		{
			if (abs(A(i,j)) > max)
			{
					max = abs(A(i,j));
					k = i;
					l = j;
			}
		}
	}
}


//Performs a rotation on a given matrix element
void jacobi(mat &A, mat &R, int k, int l, int n)
{
	double Akk = A(k,k);
	double All = A(l,l);
	double Akl = A(k,l);

	double temp;

	double t;
	double tau = (All - Akk)/(2*Akl);
	if (tau>=0) {t = 1/(tau + sqrt(1 + tau*tau));}
	else {t = 1/(tau - sqrt(1 + tau*tau));}

	double c = 1/sqrt(1 + t*t);
	double s = t*c;

	for (int i = 0; i<n; i++)
	{
		temp = A(i,k);
		A(k,i) = A(i,k) = A(i,k)*c - A(i,l)*s;
		A(l,i) = A(i,l) = A(i,l)*c + temp*s;

		temp = R(i,k);
		R(i,k) = R(i,k)*c - R(i,l)*s;
		R(i,l) = R(i,l)*c + temp*s;
	}

	double c2 = c*c;
	double s2 = s*s;

	A(k,k) = Akk*c2 - 2*Akl*c*s + All*s2;
	A(l,l) = All*c2 + 2*Akl*c*s + Akk*s2;
	A(k,l) = A(l,k) = 0;
}


//Performs Jacobis algorithm to find eigenvalues
void solveJacobi(mat A, vec &eigval, mat &eigvec, int n)
{
	int k,l;
	double max = 1;
	double eps = 1e-10;

	eigvec.zeros(n,n);
	eigvec.diag(0).fill(1.);

	while(max>eps)
	{
		max_element(A, k, l, max, n);
		jacobi(A, eigvec, k, l, n);
	}

	vec eigvalunsorted = A.diag(0);
	eigval = sort(eigvalunsorted);
	for(int i = 0; i<n; i++)
	{
		if (eigvalunsorted(i) == eigval(0))
		{
			eigvec.col(0) = eigvec.col(i);
		}
	}
}
