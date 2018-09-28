//#include <jacobi.h>
#include <iostream>
#include <math.h>
#include <armadillo>
#include <utility>
#include "jacobi.h"

using namespace arma;
using namespace std;

mat makeMatrix(vec V, double h, int n)
//----------------------------------------------------------------------------
//Creates an appropriate matrix representing the double derivative
//pluss some arbitrary potential for a given step-size
//V - potential
//h - stepsize
//n - dimmension of matrix
//----------------------------------------------------------------------------
{
	double h2 = h*h;
	mat A(n,n, fill::zeros);
	A.diag(-1).fill(-1/h2);
	A.diag(0).fill(2/h2);
	A.diag(0) += V;
	A.diag(1).fill(-1/h2);

	return A;
}

//=============================================================================
void max_element(mat &A, int &k, int &l, double &max, int n)
//----------------------------------------------------------------------------
//Find the largest off-diagonal element above the diagonal of the
//given matrix
//A - matrix to search
//k,l - indices of largest element
//max - value of largest element
//n - dimmension of matrix
//----------------------------------------------------------------------------
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


//=============================================================================
void jacobi(mat &A, mat &R, int k, int l, int n)
//----------------------------------------------------------------------------
//Does a single rotation of the matrix A to remove the element A(k,l)
//Also rotates the accompanying set of verctors R
//A - matrix to be diagonalized
//R - accompanying set of eigenvectors to be rotated
//k,l - element to be set to zero
//n - dimmension of matrix
//----------------------------------------------------------------------------
{
	double Akk = A(k,k);
	double All = A(l,l);
	double Akl = A(k,l);

	double temp;

	double t;
	double tau = (All - Akk)/(2*Akl);
	//Chooses so abs(t) is minimized
	if (tau>=0) {t = 1/(tau + sqrt(1 + tau*tau));}
	else {t = 1/(tau - sqrt(1 + tau*tau));}

	//cos and sin
	double c = 1/sqrt(1 + t*t);
	double s = t*c;

	//i may run over k,l. This will be overwritten later anyway
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


//============================================================================
int solveJacobi(mat A, vec &eigval, mat &eigvec, int n)
//----------------------------------------------------------------------------
//Diagonalize A and generate a sortet set of eigenvales and eigenvectors
//
//A - matrix to be diagonalized
//eigval - address to put eigenvalues
//eigvec - address to put eigenvectors
//n - dimmension of matrix
//showCount = true -> print total number of iterations
//----------------------------------------------------------------------------
{
	int k,l;
	double max = 1;
	double eps = 1e-10;	//Largest acceptable off-diagonal element

	mat eigvecUnsorted(n,n,fill::zeros);
	eigvecUnsorted.diag(0).fill(1.);
	eigvec.zeros(n,n);
	int count = 0;

	while(max>eps)
	{
		max_element(A, k, l, max, n);
		jacobi(A, eigvecUnsorted, k, l, n);
		count++;	//Count the number of iterations
	}

	vec eigvalUnsorted = A.diag(0);
	//The eigenvalues needs to be sorted
	eigval = sort(eigvalUnsorted);

	//Sorts the eigenvectors according to the eigenvales, since they come
	//in eigenpairs
	for(int i = 0; i<n; i++)
	{
		for(int j = 0; j<n; j++)
		{
			if (eigvalUnsorted(j) == eigval(i))
			{
				eigvec.col(i) = eigvecUnsorted.col(j);
			}
		}
	}

	return count;
}
