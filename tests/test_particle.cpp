#include "../src/particles.hpp"
#include "gtest/gtest.h"
using Eigen::RowVectorXd;
using namespace std;


// ================================================================================
// Test particle constructors
//
// ================================================================================


TEST(ParticleTests, ConstantConstructor) {
    for (int k = 1; k < 10; ++k) {
        Particle particle(k);
        EXPECT_EQ(particle.ws.size(), k);
        for (int j = 0; j < k; ++j) {
            EXPECT_DOUBLE_EQ(particle.ws(j), 1.0/k);
        }
        EXPECT_DOUBLE_EQ(particle.ws.sum(), 1.0);
    }
}


TEST(ParticleTests, VectorConstructor) {
    RowVectorXd props(3);
    props << 0.2, 0.4, 0.4;
    Particle particle(props);
    EXPECT_EQ(particle.ws, props);
}
