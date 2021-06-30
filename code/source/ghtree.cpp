
#include <ghtree.h>

#include <Dinics.h>
#include <ProgressIndicator.hpp>

using namespace std;

template<class Graph>
vector<int> residualReachable(const Graph& g, int s_or_t, bool reverse = false) {
    vector visited(size(g), false);
    visited[s_or_t] = true;
    int qB = 0; // queue begin
    vector<int> reach{s_or_t};
    while(qB<reach.size()) {
        auto v = reach[qB++];
        for(const auto& e : g[v].edges) {
            if(visited[e.to]) continue;
            if(!reverse && e.flow==e.cap) continue;
            if(reverse && -e.flow==e.cap) continue;
            visited[e.to] = true;
            reach.push_back(e.to);
        }
    }
    return reach;
}

EdgeList gomoryHuTree(const EdgeList& graph, std::vector<std::pair<int, int>> *trace) {
    Dinics::Dinics4Skip alg(graph);
    const auto& n = graph.n;

    vector parent(n, 0);
    vector flowToParent(n, numeric_limits<long long>::max());
    ProgressIndicator ind(n);
    for (int s = 1; s < n; ++s) {
        auto t = parent[s];
        if(trace) trace->emplace_back(s,t);
        alg.resetFlow();
        auto f = alg.maxFlow(s,t);
        flowToParent[s] = f;
        for(auto v : residualReachable(alg.graph(),s)) { // go through S side of cut
            // hang all siblings of s below s
            if(v!=s && parent[v]==t)
                parent[v] = s;
            // hang t below s and s below ts former parent
            if(v==parent[t]) {
                parent[s] = parent[t];
                parent[t] = s;
                flowToParent[s] = flowToParent[t];
                flowToParent[t] = f;
            }
        }
        ind.tick(s-1);
    }

    EdgeList res;
    res.n = n;
    res.hasWeights = true;
    for (int s = 1; s < n; ++s) {
        res.edges.emplace_back(s, parent[s]);
        res.weights.push_back(flowToParent[s]);
    }
    return res;
}
