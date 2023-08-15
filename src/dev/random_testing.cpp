#include "random.hpp"
#include <iostream>
#include <random>
using namespace std;

int main()
{
    RNG rng; // initialise a random number generator
    cout << "Initialise RNG with seed: " << rng.seed << endl;
    std::uniform_int_distribution<> dist(0, 10);
    for (int i = 0; i < 10; ++i) {
        cout << i << ": " << dist(rng.engine) << "\t";
    }
    cout << endl;

    RNG rng2;
    cout << "Another initialisatioin with seed: " << rng2.seed << endl;
    for (int i = 0; i < 10; ++i) {
        cout << i << ": " << dist(rng2.engine) << "\t";
    }
    cout << endl;
}