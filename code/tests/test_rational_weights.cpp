
#include <gtest/gtest.h>

#include <Dinics.h>
#include <PushRelabel.h>
#include <BoykovKolmogorov.h>

using namespace std;

TEST(rational_weights, allAlgsCanDo) {
    EdgeList g;
    g.n = 6;
    g.edges = { {0,1}, {0,2}, {1,4}, {1,3}, {2,4}, {3,5}, {4,5}, {3,2} };
    g.weights = {1.5, 2.3, 1.5, 0.25, 0.3, 1.75, 1.4, 1.35};
    g.hasWeights = true;
    g.isDirected = true;

    auto testall = [](const EdgeList& g, int s, int t, double ans) {
        ASSERT_EQ(Dinics::Dinics4Skip(g).maxFlow(s,t), ans);
        GoldbergTarjan::PushRelabel algPR(g); algPR.init(s,t);
        ASSERT_EQ(algPR.maxFlow(s,t), ans);
        ASSERT_EQ(BoykovKolmogorov::BKAlgorithm(g).maxFlow(s,t), ans);
    };

    testall(g, 0, 5, 1.65);
    g.isDirected = false;
    testall(g, 0, 5, 3.0);
    g.weights.clear();
    g.hasWeights = false;
    testall(g, 0, 5, 2.0);
}
