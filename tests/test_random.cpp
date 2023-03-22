#include <algorithm>
#include <vector>
#include "../src/random.hpp"
#include "gtest/gtest.h"


// ================================================================================
// Test random number generator seed creation
//
// ================================================================================


TEST(RandomNumberGeneratorTests, CheckSeedsNotDuplicated) {
    int N = 1000;
    std::vector<int> random_seeds(N);
    for (int i = 0; i < N; ++i) {
        RNG rng;                                // Should initialise and store seed
        random_seeds[i] = rng.seed;
    }
    std::sort(random_seeds.begin(), random_seeds.end());
    auto found = std::adjacent_find(random_seeds.begin(), random_seeds.end());
    EXPECT_EQ(found, random_seeds.end());
}

