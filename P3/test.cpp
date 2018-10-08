#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <vector>
using namespace std;
using namespace arma;

int main()
{
    cube pos(3,4,8,fill::zeros);
    cout << pos.slice(0) << endl;
    return 0;
}
