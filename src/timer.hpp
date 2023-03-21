#pragma once
#include <chrono>
using namespace std;

// TimeUnit could be:
// chrono::nanoseconds
// chrono::microseconds
// chrono::miliseconds
// chrono::seconds
// &c.
template<typename TimeUnit>
class Timer
{
private:
    using Clock = chrono::steady_clock;
    chrono::time_point<Clock> start_time;

public:
    // Constructors
    Timer() {};

    // Methods
    void start() {
        start_time = Clock::now();
    }
    double elapsed() {
        return chrono::duration_cast<TimeUnit>(Clock::now() - start_time).count();
    }
};