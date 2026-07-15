#include <vector>
#include <iostream>
using namespace std;

struct Params{
    bool biseparable;
    unsigned int maxmodulo;
    long double GCweight;
};

vector<string> generate_seeds(string structure, unsigned int target_nb_sols, Params* p);
/*
    Input:
        * structure, a secondary structure as a well-parenthesized word
        * target_nb_sols, the target number of solutions
        * p, a list of parameters including:
            * biseparable, boolean, 1 if unpaired positions feature both A and C
            * maxmodulo, the maximal modulus to be considered for m/mA/mC
            * GCweight, the weight of the GC pairs, default is 0.0
            
    Output:
        * A list of designs of size at most target_nb_sols
*/