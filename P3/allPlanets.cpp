#include <iostream>
#include <iomanip>
#include <armadillo>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include "solver.h"
#include "planet.h"
using namespace std;
using namespace arma;

double scale = 4*M_PI*M_PI;

int main()
{
    ifstream myfile;
    myfile.open("initial_conditions.txt");
    string line;

    //vector<Planet> solarSystem;

    if(myfile.is_open())
    {
        //myfile.ignore(256, '\n');   //Ignores first line
        while(getline(myfile,line))
        {
            istringstream buf(line);
            istream_iterator<string> beg(buf), end;
            vector<string> tokens(beg, end);

            //vec posIC = {stof(tokens[1]),stof(tokens[2]),stof(tokens[3])};
            //vec velIC = {stof(tokens[4]),stof(tokens[5]),stof(tokens[6])};
            //velIC *= 365.25;
            //double mass = stof(tokens[7]);
            //Planet planet(posIC, velIC, mass);

            //solarSystem.push_back(planet);
        }
    }
    return 0;
}
