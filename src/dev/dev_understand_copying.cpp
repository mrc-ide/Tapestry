#include <iostream>
#include "libs/eigen-3.4.0/Dense"
#include "timer.hpp"
using Eigen::MatrixXd;
using namespace std;


class NaiveWay
{
private:
    MatrixXd array;

public:
    NaiveWay(int K, int N)
    {
        array = MatrixXd::Constant(K, N, 0.0);
        double z = 0.0;
        for (int i = 0; i < K; ++i) {
            for (int j =0 ; j < N; ++j) {
                array(i, j) = z;
            }
        }
    }
};

// This is probably fastest, across different matrix sizes
class MemberInitList
{
private:
    MatrixXd array;

public:
    MemberInitList(int K, int N)
        : array(MatrixXd::Constant(K, N, 0.0))
    {
        double z = 0.0;
        for (int i = 0; i < K; ++i) {
            for (int j =0 ; j < N; ++j) {
                array(i, j) = z;
            }
        }
    }
};

// Probably slowest
class StaticMethod
{
private:
    const MatrixXd array;
    MatrixXd static create_array(int K, int N)
    {
        MatrixXd eg_array = MatrixXd::Constant(K, N, 0.0);
        double z = 0.0;
        for (int i = 0; i < K; ++i) {
            for (int j =0 ; j < N; ++j) {
                eg_array(i, j) = z;
            }
        }
        return eg_array;
    }
public:
    StaticMethod(int K, int N)
        : array(create_array(K, N))
    {};
};


int main(int argc, char* argv[])
{
    // Fixed
    int n_iters = 10000;

    // User
    if (argc < 3) {
        cout << "Incorrect usage: ./aout <int> <int>" << endl;
        exit(1);
    }
    int K = atoi(argv[1]);
    int N = atoi(argv[2]);
    
    // TESTING
    // ---------------------------------------------------------------------------
    cout << "Naive Way";
    Timer<chrono::microseconds> timer;
    timer.start();
    for (int i = 0; i < n_iters; ++i) {
        NaiveWay naive_way(K, N);
    };
    cout << "  Time elapsed (us):" << timer.elapsed() << endl;

    // ---------------------------------------------------------------------------
    cout << "Member Init List";
    timer.start();
    for (int i = 0; i < n_iters; ++i) {
        MemberInitList member_init(K, N);
    };
    cout << "  Time elapsed (us):" << timer.elapsed() << endl;

    // ---------------------------------------------------------------------------
    cout << "Static Approach";
    timer.start();
    for (int i = 0; i < n_iters; ++i) {
        StaticMethod static_method(K, N);
    };
    cout << "  Time elapsed (us):" << timer.elapsed() << endl;
}