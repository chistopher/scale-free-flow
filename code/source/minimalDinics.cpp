
#include <minimalDinics.h>

#include <queue>
#include <cassert>

#define rep(a, b)   for(int a = 0; a < (b); ++a)

using namespace std;
using namespace MinimalDinics;
const int INF = 1e9;

// DFS for normal distances
int dfs(int v, int aug, int t,
        vector<int>& next,
        vector<vector<edge>>& adj,
        const vector<int>& dist)
{
    if (v == t) return aug;
    for (int& i = next[v]; i<adj[v].size(); ++i) {
        auto& e = adj[v][i];
        if (e.flow == e.cap) continue;
        if (dist[e.to] != dist[v] + 1) continue;
        int pushed = dfs(e.to, min(aug, e.cap - e.flow), t, next, adj, dist);
        if (pushed == 0) continue;
        e.flow += pushed;
        adj[e.to][e.rev].flow -= pushed;
        return pushed;
    }
    return 0;
}

// computes all distances up to the layer of t
vector<int> bfs(int s, int t, const vector<vector<edge>>& adj) {
    vector<int> dist(adj.size(),INF);
    dist[s] = 0;
    queue<int> q{{s}}; // TODO use vector
    while(!q.empty()) {
        auto v = q.front(); q.pop();
        if(dist[v] == dist[t]) break;
        for(auto& e : adj[v]) {
            if(dist[e.to] == INF && e.cap > e.flow) {
                dist[e.to] = dist[v] + 1;
                q.push(e.to);
            }
        }
    }
    return dist;
}

long long MinimalDinics::dinics(int s, int t, std::vector<std::vector<edge>>& adj) {
    auto flow =0ll;
    while(true) {
        auto dist = bfs(s,t,adj);
        if(dist[t]==INF) break;
        vector<int> next(size(adj),0);
        int aug;
        while((aug = dfs(s,INF,t,next,adj,dist))) flow += aug;
    }
    return flow;
};

vector<vector<edge>> MinimalDinics::buildNetwork(const vector<pair<int,int>>& edges, int n) {
    vector<vector<edge>> adj(n);
    for(auto [a,b] : edges) {
        assert(0<=a && a<n);
        assert(0<=b && b<n);
        adj[a].push_back(edge{a,b,0,1,adj[b].size()});
        adj[b].push_back(edge{b,a,0,1,adj[a].size()-1});
    }
    return adj;
}

vector<vector<edge>> MinimalDinics::buildNetwork(const vector<array<int,3>>& edges, int n) {
    vector<vector<edge>> adj(n);
    for(auto [a,b,w] : edges) {
        assert(0<=a && a<n);
        assert(0<=b && b<n);
        adj[a].push_back(edge{a,b,0,w,adj[b].size()});
        adj[b].push_back(edge{b,a,0,w,adj[a].size()-1});
    }
    return adj;
}

vector<int> MinimalDinics::cutSide(int s, int t, const vector<vector<edge>>& adj, bool sideS) {
    vector<bool> visited(size(adj),false);
    visited[sideS ? s : t] = true;
    int qB = 0; // queue begin
    vector<int> reach{sideS ? s : t};
    while(qB<reach.size()) {
        auto v = reach[qB++];
        for(auto& e : adj[v]) {
            if(visited[e.to]) continue;
            if(sideS && e.flow==e.cap) continue;
            if(!sideS && adj[e.to][e.rev].flow==adj[e.to][e.rev].cap) continue;
            visited[e.to] = true;
            reach.push_back(e.to);
        }
    }
    return reach;
}
