#include <func.h>
using namespace arma;
using namespace std;

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
void Jacobi(mat &A, int k, int l, int n)
{
	double Akk = A(k,k);
	double All = A(l,l);
	double Akl = A(k,l);

	double t;
	double tau = (All - Akk)/(2*Akl);
	if (tau>=0) {t = 1/(tau + sqrt(1 + tau*tau));}
	else {t = 1/(tau - sqrt(1 + tau*tau));}

	double c = 1/sqrt(1 + t*t);
	double s = t*c;

	for (int i = 0; i<n; i++)
	{
		double temp = A(i,k);
		A(k,i) = A(i,k) = A(i,k)*c - A(i,l)*s;
		A(l,i) = A(i,l) = A(i,l)*c + temp*s;
	}

	double c2 = c*c;
	double s2 = s*s;

	A(k,k) = Akk*c2 - 2*Akl*c*s + All*s2;
	A(l,l) = All*c2 + 2*Akl*c*s + Akk*s2;
	A(k,l) = A(l,k) = 0;
}


//Performs Jacobis algorithm to find eigenvalues
vec SolveJacobi(double h, vec &x, vec V(double, vec), double w, int n)
{
	double h2 = h*h;
	mat A(n,n, fill::zeros);
	A.diag(-1).fill(-1/h2);
	A.diag(0).fill(2/h2);
	A.diag(0) += V(w,x);
	A.diag(1).fill(-1/h2);

	int k,l;
	double max = 1;
	double eps = 1e-10;

	while(max>eps)
	{
		max_element(A, k, l, max, n);
		Jacobi(A, k, l, n);
	}

	vec eig = A.diag(0);
	eig = sort(eig);
	return eig;
}

mat SolveArma(int n, double h, vec &x, vec V(double, vec), double w)
{
	double h2 = h*h;
	mat A(n,n, fill::zeros);
	A.diag(-1).fill(-1/h2);
	A.diag(0).fill(2/h2);
	A.diag(0) += V(w,x);
	A.diag(1).fill(-1/h2);

	vec eigval;
	mat eigvec;

	eig_sym(eigval, eigvec, A);
	return eigvec;
}
