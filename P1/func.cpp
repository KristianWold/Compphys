#include <iostream>
#include <chrono>
#include <math.h>

double **matrix(int, int);
double *vector(int);

double **matrix(int row, int col)
{
  double ** A;
  A = new double*[row];
  for (int i = 0; i < row; i++)
  {
    A[i] = new double[col];
  }
  return (double **)A;
}

double *vector(int row)
{
  double * A;
  A = new double[row];
  return (double *)A;
}

void **free_matrix(int n, double** A)
{
  for (int i = 0; i<n; i++)
  {
    delete[] A[i];
  }
  delete[] A;
}

void *free_vector(double* A)
{
  delete[] A;
}

double source(double x){return 100*exp(-10*x);}

double solution(double x){return 1 - (1-exp(-10))*x - exp(-10*x);}


//General method. Solves the set of linear equations with a forward and backward substitution
double *solveB(int n, double source(double))
{
  double h = 1/double(n-1);

  double *a = vector(n); double *b = vector(n); double *c = vector(n);
  for(int i = 0; i<n; i++)
  {
    a[i] = -1.;
    b[i] = 2.;
    c[i] = -1.;
  }
  //source term
  double* f = vector(n);
  for(int i = 0; i<n; i++)
  {
    f[i] = pow(h,2)*source(i*h);
  }

  auto start = std::chrono::high_resolution_clock::now(); //start clock
  //forward sub
  double temp;
  for(int i=2; i<n-1; i++)
  {
    temp = a[i-1]/b[i-1];
    b[i] -= temp*c[i-1];
    f[i] -= temp*f[i-1];
  }

  //backward sub
  double *v = vector(n);
  v[0] = 0; v[n-1] = 0; //boundary conditions
  for(int i = n-2; i>0; i--)
  {
    v[i] = (f[i] - c[i]*v[i+1])/b[i];
  }

  auto finish = std::chrono::high_resolution_clock::now(); //stop clock
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time:" << elapsed.count() << std::endl;

  return (double *)v;
  free_vector(a); free_vector(b); free_vector(c);
  free_vector(v); free_vector(f);
}

//Optimized solver
double *solveC(int n, double source(double))
{
  double h = 1/double(n-1);

  double *b = vector(n);
  for(int i = 0; i<n; i++)
  {
    b[i] = 2.;
  }

  //source term
  double* f = vector(n);
  for(int i = 0; i<n; i++)
  {
    f[i] = pow(h,2)*source(i*h);
  }

  auto start = std::chrono::high_resolution_clock::now(); //start clock
  //forward sub
  for(int i=2; i<n-1; i++)
  {
    b[i] -= 1/b[i-1];
    f[i] += f[i-1]/b[i-1];
  }

  //backward sub
  double *v = vector(n);
  v[0] = 0; v[n-1] = 0; //boundary condition
  for(int i = n-2; i>0; i--)
  {
    v[i] = (f[i] + v[i+1])/b[i];
  }

  auto finish = std::chrono::high_resolution_clock::now(); //stop clock
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time:" << elapsed.count() << std::endl;

  return (double *)v;
  free_vector(b); free_vector(v); free_vector(f);
}
