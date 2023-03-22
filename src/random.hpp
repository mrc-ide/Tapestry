#pragma once
#include <random>


/*
* Simple class to wrap the generation of random numbers and store seed
* TODO:
* - Important to confirm this is best-practice
*/
struct RNG
{
    const int seed;
    std::mt19937 engine;

    RNG()
        : seed([]() -> int { std::random_device rd; return rd();}()),
        engine(seed)
    {}
};
