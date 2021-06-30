
#include <gtest/gtest.h>

#include "global_network_data.hpp"

#include <BoykovKolmogorov.h>

// ACTUAL TESTS
TEST(BoykovKolmogorovTest, randomPairs) {
    BoykovKolmogorov::BKAlgorithm alg(testGraph);
    for (auto [s, t, ans] : pairs) {
        EXPECT_EQ(alg.maxFlow(s, t), ans);
    }
}
