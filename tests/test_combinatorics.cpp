#include <stdexcept>
#include <numeric>
#include "../src/combinatorics.hpp"
#include "../src/constants.hpp"
#include "gtest/gtest.h"
using namespace std;


// ================================================================================
// Test creation of powersets
//
// ================================================================================


TEST(PowersetTests, InputExceptions) {
    EXPECT_THROW(create_powerset(-1), std::invalid_argument);
    EXPECT_THROW(create_powerset(0), std::invalid_argument);
    EXPECT_NO_THROW(create_powerset(1));
}

TEST(PowersetTests, MatrixSize) {
    int MAX_K = 10;
    for (int k = 1; k < MAX_K; ++k) {
        int N = std::pow(2, k);
        MatrixXi powerset = create_powerset(k);
        EXPECT_EQ(powerset.rows(), N);
        EXPECT_EQ(powerset.cols(), k);
    }
}

TEST(PowersetTests, MatrixSum) {
    int MAX_K = 10;
    for (int k = 1; k < MAX_K; ++k) {
        MatrixXi powerset = create_powerset(k);
        int expected_sum = k * std::pow(2, k - 1); // each element is include half of time
        EXPECT_EQ(powerset.sum(), expected_sum);
    }
}

int calc_factorial(const int& value) {
        int result = 1;
        for (int i = value; i > 0; --i) {
            result *= i;
        }
        return result;
}

TEST(PowersetTests, RowSums) {
    int MAX_K = 10;
    for (int k = 1; k < MAX_K; ++k) {
        MatrixXi powerset = create_powerset(k);
        Eigen::VectorXi row_sums = powerset.rowwise().sum();
        for (int j = 0; j <= k; ++j) {
            int expected_with_j = calc_factorial(k) / (calc_factorial(k - j) * calc_factorial(j));
            int found_with_j = 0;
            for (const int& r : row_sums) {
                if (r == j) {
                    ++found_with_j;
                }
            }
            EXPECT_EQ(expected_with_j, found_with_j);
        }
    }
}


// ================================================================================
// Test creation of all possible partitions
//
// ================================================================================


TEST(PartitionTests, NumberOfPartitions) {
    int MAX_K = 10;
    for (int k = 1; k < MAX_K; ++k) {
        vector<int> collection(k);
        std::iota(collection.begin(), collection.end(), 0);
        vector<vector<vector<int>>> partitions = create_all_partitions(collection);
        EXPECT_EQ(partitions.size(), BELL_NUMBERS[k]);
    }
}

