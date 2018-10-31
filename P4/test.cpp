#include <iostream>
#include <map>
#include <armadillo>
#include <random>
using namespace std;
using namespace arma;

class Ost
{
public:
    int L;
    int test;
    uniform_int_distribution<int> rand_coord(0,1);
    Ost(int L_, mt19937 &engine)
    {
        //rand_coord = uniform_int_distribution<int>(0,L_-1);
        test = rand_coord(engine);
    }
};

int main(int argc, char const *argv[])
{
    mt19937 engine(1);
    Ost parmesan(10, engine);

    return 0;
}
