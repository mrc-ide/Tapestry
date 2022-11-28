#include <vector>
#include <numeric>
#include <cmath>
#include <testthat.h>
#include "combinatorics.h"
using namespace std;


/**
 * Testing creation of powersets
 * 
 * To include:
 * - Bad inputs
 * - Sums per row / column
 * - Something to ensure orientation is OK
 * 
*/
context("Binary powerset matrix")
{
    // Create range of COI values over which to test
    int MAX_K = 10;
    vector<int> Ks(MAX_K - 1);
    std::iota(Ks.begin(), Ks.end(), 1);

    test_that("Test size of powerset matrix")
    {

        for (int k : Ks)
        {
            int N = std::pow(2, k);
            vector<vector<int>> powerset = create_powerset(k);

            expect_true(powerset.size() == N);
            expect_true(powerset[0].size() == k);
        }
    }

    test_that("Test contents sum to expected")
    {

        for (int k : Ks)
        {
            vector<vector<int>> powerset = create_powerset(k);

            int sum = 0;
            for (int i = 0; i < powerset.size(); ++i)
            {
                for (int j = 0; j < powerset[0].size(); ++j)
                {
                    sum += powerset[i][j];
                };
            };

            // Each element should be included HALF of the time
            int expect_sum = k * std::pow(2, k - 1);
            expect_true(sum == expect_sum);
        }
    }
}


/**
 * Testing creation of all possible partitions
 * 
 * To include:
 * - Can I ensure the sizes are right?
 * - Can I ensure I always find certain specific partitions?
 * 
*/
context("Test all possible partitions")
{
    // Create range of COI values over which to test
    int MAX_K = 10;
    vector<int> Ks(MAX_K - 1);
    std::iota(Ks.begin(), Ks.end(), 1);

    // One to ten; from https://oeis.org/A000110
    vector<int> expected_bell_numbers = {
        1,
        1,
        2,
        5,
        15,
        52,
        203,
        877,
        4140,
        21147,
        115975};

    test_that("Test correct number of partitions, i.e. Bell(k)")
    {
        for (int k : Ks)
        {
            // Create a collection of objects
            vector<int> collection(k);
            std::iota(collection.begin(), collection.end(), 0);

            // Create all possible partitions of them
            vector<vector<vector<int>>> partitions; 
            partitions = create_all_partitions(collection);

            // Check there are a Bell number of partitions
            expect_true(partitions.size() == expected_bell_numbers[k]);
        }
    }
}