#pragma once
#include <chrono>
#include <ctime>
#include <sstream>
#include <string>
using namespace std;

// TimeUnit could be:
// chrono::nanoseconds
// chrono::microseconds
// chrono::milliseconds
// chrono::seconds
// &c.
class Timer
{
private:
    using Clock = chrono::steady_clock;
    chrono::time_point<Clock> start_time;

public:
    // Constructor
    Timer()
        : start_time(Clock::now())
        {};


    // Reset the timer
    void reset() {
        start_time = Clock::now();
    }

    // Return time elapsed
    template<typename TimeUnit>
    double elapsed() {
        return chrono::duration_cast<TimeUnit>(Clock::now() - start_time).count();
    }

    // Print start time
    // TODO: can't figure this out, will take like a day of screwing around with chrono
    // string get_current_time()
    // {
    //     //chrono::time_point<Clock> now =  Clock::now();
    //     time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    //     return std::to_string(ctime(&now));
    // }
};