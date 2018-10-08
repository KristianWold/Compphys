#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <vector>
using namespace std;
using namespace arma;

int main()
{
    cube pos(30,40,8,fill::zeros);
    cout << &pos.slice(0) << endl;
    cout << &pos.slice(1) << endl;
    cout << &pos.slice(2) << endl;
    return 0;
}
