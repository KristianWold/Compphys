#include<cmath>
#include<random>
#include<iostream>
#include<iomanip>
#include<chrono>
#include<armadillo>

using namespace std;
using namespace chrono;
using namespace arma;

typedef double (*myfunc)(mat);

class VMC
{
public:
    myfunc acceptAmp;
    VMC(myfunc newAcceptAmp)
    {
        acceptAmp = newAcceptAmp;
    }
};

double acceptAmp(mat r)
{
    return r(0,0);
}

int main(int argc, char const *argv[]) {
    VMC solver(&acceptAmp);
    mat r = ones(3,3);
    cout << solver.acceptAmp(r);
    return 0;
}
