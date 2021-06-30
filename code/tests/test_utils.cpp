
#include <gtest/gtest.h>

#include <random>

#include <utils.h>

using namespace std;

mt19937 gen(1337);

TEST(utils, CConeCircle) {
    EdgeList g1;
    g1.n = 100;
    for(int i=0; i<100; ++i) g1.edges.emplace_back(i, (i+1)%g1.n);
    auto g2 = largestCC(g1);
    ASSERT_EQ(g1.edges,g2.edges);
}

TEST(utils, CCgaplessRemapping) {
    EdgeList g;
    g.n = 100;
    uniform_int_distribution<int> dist(0,g.n-1);
    for(int i=0; i<g.n-2; ++i)
        g.edges.emplace_back(dist(gen),dist(gen));
    auto g2 = largestCC(g);
    EXPECT_NE(g.n,g2.n);
    EXPECT_NE(g.edges,g2.edges);
    vector<bool> used(g2.n,false);
    for(auto e : g2.edges) {
        for(auto v : {e.first, e.second}) {
            ASSERT_GE(v,0);
            ASSERT_LT(v,g.n);
            used[v] = true;
        }
    }
    for(int v=0; v<g2.n; ++v)
        ASSERT_TRUE(used[v]);
}

TEST(utils, degreeSanity) {
    int n = 100;
    EdgeList g;
    g.n = n;
    uniform_int_distribution dist(0,n-1);
    for(int i=0; i<10*n; ++i)
        g.edges.emplace_back(dist(gen),dist(gen));
    auto deg = degrees(g);
    for(auto d : deg) {
        ASSERT_GE(d,0);
        ASSERT_LE(d,size(g.edges));
    }
    ASSERT_EQ(2*size(g.edges),accumulate(begin(deg),end(deg),0ll));
}
