#include<armadillo>

using namespace arma;
/*
int main()
{
    mat A = ones(3,3);
    vec a = sum(A%A,1);
    cout << sum(a) << endl;
    return 0;
}
*/

int main()
{
    vec A = -1*ones(3);
    cout << norm(A,2) << endl;
    return 0;
}
