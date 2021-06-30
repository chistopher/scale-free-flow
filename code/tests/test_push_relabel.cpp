
#include <gtest/gtest.h>

#include "global_network_data.hpp"

#include <PushRelabel.h>

// ACTUAL TESTS
TEST(PushRelabelTest, randomPairs) {
    GoldbergTarjan::PushRelabel alg(testGraph);
    for (auto [s, t, ans] : pairs) {
        alg.init(s,t);
        EXPECT_EQ(alg.maxFlow(s, t), ans);
        alg.deinit();
    }
}
