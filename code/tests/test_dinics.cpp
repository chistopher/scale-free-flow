
#include <gtest/gtest.h>

#include "global_network_data.hpp"

#include <bitset>

#include <Dinics.h>


template<class DinicsType>
class DinicsTest : public ::testing::Test { };

template<unsigned int d, class ... ALGs>
struct Combinations {
    static_assert(d != 0);
    constexpr static auto b = std::bitset<5>(d);
    using type = typename Combinations<d - 1, Dinics::Dinics<b[0], b[1], b[2], b[3], b[4]>, ALGs...>::type;
};
template<class ... ALGs>
struct Combinations<0, ALGs...> {
    using type = ::testing::Types<Dinics::Dinics<0, 0, 0, 0, 0>, ALGs...>;
};

using AllAlgs = Combinations<31>::type;
//using AllAlgs = testing::Types<Dinics2<true, true, true, true, true>>;
TYPED_TEST_SUITE(DinicsTest, AllAlgs);

// ACTUAL TESTS
TYPED_TEST(DinicsTest, randomPairs) {
    TypeParam alg(testGraph);
    for (auto [s, t, ans] : pairs) {
        EXPECT_EQ(alg.maxFlow(s, t), ans);
        alg.resetFlow();
    }
}

constexpr std::array statResults{
    "3,0,3,440,525,7080,2760,0,0,0,0,0,0,0,100,132,0,0,1753,795,0,0"  "0,3,7,989,0,9893,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0", // uni nobreak noskip
    "3,0,3,440,525,7080,2760,0,0,0,0,0,0,0,100,0,0,0,1753,0,0,0"      "0,3,7,989,0,9893,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0", // uni nobreak skip
    "3,0,3,440,53,7080,280,0,0,0,0,0,0,0,100,52,0,0,1753,277,0,0"     "0,3,7,989,0,9893,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0", // uni break noskip
    "3,0,3,440,53,7080,280,0,0,0,0,0,0,0,100,0,0,0,1753,0,0,0"        "0,3,7,989,0,9893,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0", // uni break skip

    "3,0,1,1,5,5,545,2,4,193,284,4455,4,541,1,5,4,3,5,106,102,264"    "0,3,0,0,1,0,5,1,1,0,3,0,0,0,0,0,0,0,0,0,0,0", // bi noskip
    "3,0,1,1,5,5,545,2,4,193,284,4455,4,541,1,4,4,3,5,102,102,264"    "0,3,0,0,1,0,5,1,1,0,3,0,0,0,0,0,0,0,0,0,0,0", // bi skip
};

TYPED_TEST(DinicsTest, hasCorrectSearchSpace) {
    TypeParam alg(testGraph);
    auto [s,t,w] = pairs.front();
    auto f = alg.maxFlow(s,t,true);
    EXPECT_EQ(f,w);
    std::stringstream ss;
    for(auto& round : alg.stats()) ss << round;
    auto resID = 4*TypeParam::isBI + 2*(!TypeParam::isBI && TypeParam::isADJREV) + TypeParam::isSKIP;
    ASSERT_EQ(statResults[resID], ss.str());
};


