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
    mat A = ones(3,2);
    cout << A.col(0) << endl;
    return 0;
}
