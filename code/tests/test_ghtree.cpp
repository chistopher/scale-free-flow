
#include <gtest/gtest.h>

#include <ghtree.h> // being tested
#include <Dinics.h>
#include <utils.h>

using namespace std;

#ifndef NDEBUG
const EdgeList smallGraph{ {
        {0,1}, {0,2}, {0,3},
        {3,4}, {2,3},
        {0,5}, {1,5},
        {5,6}, {5,7}, {5,8},
        {6,7},
}, {}, false, false, 9};

/*
 *      1   2 - 3 - 4
 *      | \ | /
 *      5 - 0
 *    / | \
 *  6 - 7   8
 */
#else
#include "global_network_data.hpp"
const auto& smallGraph = testGraph;
#endif


TEST(ghtree, correctNumberOfEdges) {
    auto tree = gomoryHuTree(smallGraph);
    ASSERT_EQ(tree.n, smallGraph.n);
    ASSERT_EQ(size(tree.edges), smallGraph.n-1);
}

TEST(ghtree, writesTrace) {
    vector<pair<int,int>> trace;
    auto tree = gomoryHuTree(smallGraph, &trace);
    ASSERT_EQ(size(trace), size(tree.edges));
}

auto asAdjacency(const EdgeList& tree) {
    vector<vector<pair<int,double>>> adj(tree.n);
    for(int i=0; i<size(tree.edges); ++i) {
        auto [a,b] = tree.edges[i];
        auto w = tree.weights[i];
        adj[a].emplace_back(b,w);
        adj[b].emplace_back(a,w);
    }
    return adj;
}

// returns weighted edge in ghtree
auto minOnPath(const vector<vector<pair<int,double>>>& adj, int s, int t) {
    tuple<int,int,double> res{0,0,1e9};
    auto dfs = [&](int cur, int par, int target, auto&& recur) {
        if(cur==target) return true;
        for(auto [nei, w] : adj[cur]) {
            if(nei==par) continue;
            auto found = recur(nei,cur,target,recur);
            if(!found) continue;
            if(get<2>(res)>w) res = {cur,nei,w};
            return true;
        }
        return false;
    };
    dfs(s,-1,t,dfs);
    return res;
}

TEST(ghtree, flowEquivalent) {
    const auto n = smallGraph.n;
    auto tree = gomoryHuTree(smallGraph);
    auto adj = asAdjacency(tree);
    Dinics::Dinics5OPT alg(smallGraph);
    for (int s = 0; s < n; ++s) {
        for (int t = s + 1; t < n; ++t) {
            auto edge = minOnPath(adj, s, t);
            alg.resetFlow();
            ASSERT_EQ(get<2>(edge), alg.maxFlow(s,t));
        }
    }
}

double cutValue(const vector<bool>& inS, const EdgeList& graph) {
    auto ans = 0.0;
    for(int i = 0; i<size(graph.edges); ++i) {
        auto [u,v] = graph.edges[i];
        auto w = graph.hasWeights ? graph.weights[i] : 1.0;
        if(inS[u] == inS[v]) continue;
        if(!graph.isDirected || (inS[u] && !inS[v]))
            ans += w;
    }
    return ans;
}

void nodesBehind(int v, int par, const vector<vector<pair<int,double>>>& adj, vector<bool>& behind) {
    behind[v] = true;
    for(auto [nei,_] : adj[v])
        if(nei!=par) nodesBehind(nei,v,adj,behind);
}

TEST(ghtree, cutEquivalent) {
    const auto n = smallGraph.n;
    auto tree = gomoryHuTree(smallGraph);
    auto adj = asAdjacency(tree);
    for (int s = 0; s < n; ++s) {
        for (int t = s + 1; t < n; ++t) {
            auto [u,v,w] = minOnPath(adj, s, t);
            vector<bool> inS(n, false);
            nodesBehind(u,v,adj,inS);
            ASSERT_TRUE(inS[s]);
            ASSERT_TRUE(inS[u]);
            ASSERT_FALSE(inS[v]);
            ASSERT_FALSE(inS[t]);
            ASSERT_EQ(cutValue(inS,smallGraph), w);
        }
    }
}
