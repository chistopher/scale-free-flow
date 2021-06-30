#include <bits/stdc++.h>

#include <omp.h>
#include <girgs/Generator.h>

#define rep(a, b)   for(int a = 0; a < (b); ++a)
#define all(a)      (a).begin(),(a).end()
//#define endl        '\n'

using namespace std;
int INF = 1e9;


void make_pdf(
        const std::vector<std::vector<double>> &pos,
        const std::vector<std::pair<int, int>> &graph,
        std::string file,
        map<int,string> nodeColors = {},
        map<int,string> nodeLabels = {},
        map<pair<int,int>,string> edgeColors = {}) {

    std::ofstream f{file+".dot"};
    //f << "graph girg {\n\toverlap=scale;\n\n";
    f << "graph girg {\n\tscale=" << 2000 << ";\n\n";
    f << std::fixed;
    for (int i = 0; i < pos.size(); ++i) {
        f << '\t' << i << " [pos=\"" << pos[i][0] << ',' << pos[i][1] << "\"";
        if(nodeColors.count(i) > 0)
            f << " style =\"filled\" fillcolor=\"" << nodeColors[i] << "\"";
        if(nodeLabels.count(i) > 0)
            f << " label=\"" << nodeLabels[i] << "\"";
        f << "];\n";
    }
    f << '\n';
    for (auto e : graph) {
        auto [a,b] = minmax(e.first,e.second);
        f << '\t' << a << "\t-- " << b;
        if(edgeColors.count({a,b}) > 0)
            f << "\t[color=\"" << edgeColors[{a,b}] << "\"]";
        f << ";\n";
    }
    f << "}\n";
    f.flush();

    // create pdf
    string command = "neato -n -Tpdf " + file + ".dot -o " + file + ".pdf";
    system(command.c_str());
}




struct edge {
    int from, to;
    int flow, cap;
    size_t rev;
};

int dfs(int v, int aug, int t, vector<int>& dist, vector<int>& next, vector<vector<edge>>& adj) {
    if (v == t) return aug;
    for (int& i = next[v]; i<adj[v].size(); ++i) {
        auto e = &adj[v][i];
        if (e->flow == e->cap) continue;
        if (dist[e->to] != dist[v] + 1) continue;
        int pushed = dfs(e->to, min(aug, e->cap - e->flow), t, dist, next, adj);
        if (pushed == 0) continue;
        e->flow += pushed;
        adj[e->to][e->rev].flow -= pushed;
        return pushed;
    }
    return 0;
}

vector<vector<double>> pos;
vector<pair<int, int>> edges;

vector<int> bfs_dist(int s, int t, int n, vector<vector<edge>>& adj) {
    set<int> space;
    set<pair<int,int>> edge_space;

    vector<int> dist(n,INF);
    dist[s] = 0;
    queue<int> q{{s}};
    while(!q.empty()) {
        auto v = q.front(); q.pop();
        if(dist[v] == dist[t]) break;
        space.insert(v);
        for(auto& e : adj[v]) {
            edge_space.insert(minmax(v,e.to));
            if(dist[e.to] == INF && e.cap > e.flow) {
                dist[e.to] = dist[v] + 1;
                q.push(e.to);
            }
        }
    }

    cout << "bfs\t" << space.size() << '\t' << edge_space.size() << endl;
    return dist;
}

vector<int> bfs_dist_bi(int s, int t, int n, vector<vector<edge>>& adj) {

    map<int,string> total_space{{s,"red"},{t,"blue"}};

    vector<int> dist_s(n,INF), dist_t(n,INF);
    deque<int> qs{{s}}, qt{{t}};
    dist_s[s] = 0; dist_t[t] = 0;
    int shortest = INF;
    auto process = [&](bool reverse) {
        auto& dist = reverse ? dist_t : dist_s;
        auto& dist_other = !reverse ? dist_t : dist_s;
        auto& q = reverse ? qt : qs;
        if(q.empty()) return;
        auto v = q.front(); q.pop_front();
        if((reverse&&v==s)||(!reverse&&v==t)) return;

        if(!total_space.count(v)) total_space[v] = reverse ? "blue" : "red";
        map<int,string> space = total_space;
        space[v] = reverse ? "blue4" : "red4";
        map<pair<int,int>,string> edge_space;

        for(auto& e : adj[v]) {
            if(!total_space.count(e.to)) space[e.to] = "green4";
            edge_space[{minmax(v,e.to)}] = "yellow";
            auto full = (reverse ? adj[e.to][e.rev].flow == adj[e.to][e.rev].cap : e.flow==e.cap);
            if(full || dist[e.to] <= dist[v]+1) continue;
            edge_space[{minmax(v,e.to)}] = "green";
            dist[e.to] = dist[v]+1;
            if(dist_other[e.to] < INF) shortest = min(shortest, dist[e.to]+dist_other[e.to]);
            if(shortest == INF || (shortest!=INF && dist[e.to]+dist_other[e.to] == shortest)) {
                space[e.to] = "green";
                q.push_back(e.to);
            }
        }

        for(auto u : q) if(!space.count(u)) space[u] = "gray";
        make_pdf(pos,edges,"view",space,{},edge_space);
        auto a = 0;
    };

    while(qs.size() || qt.size()) {
        if(qs.size()) process(false);
        if(shortest==INF&&qs.empty()) break;
        if(qt.size()) process(true);
    }

    make_pdf(pos,edges,"view",total_space);

    rep(i,n) if(dist_s[i]+dist_t[i]!=shortest) dist_s[i] = INF;

    return dist_s;
}

vector<int> bfs_dist_bi2(int s, int t, int n, vector<vector<edge>>& adj) {

    set<int> space;
    set<pair<int,int>> edge_space;

    vector<int> dist_s(n,INF), dist_t(n,INF);
    dist_s[s] = 0; dist_t[t] = 0;
    deque<int> qs{{s}}, qt{{t}};
    bool found = false;

    // DEBUG
    auto print = [&]{
        int sLayer = 0, tLayer = 0;
        rep(i,n) if(dist_s[i]!=INF) sLayer = max(sLayer,dist_s[i]);
        rep(i,n) if(dist_t[i]!=INF) tLayer = max(tLayer,dist_t[i]);
        map<int,string> color{{s,"red"},{t,"blue"}};
        rep(i,n) {
            auto [ds,dt] = pair(dist_s[i],dist_t[i]);
            if(ds==INF&&dt==INF) continue;
            if(ds<sLayer) color[i] = "red";
            else if(dt<tLayer) color[i] = "blue";
            else if(ds==sLayer&&dt==tLayer) color[i] = "yellow";
            else if(ds==sLayer) color[i] = "red4";
            else if(dt==tLayer) color[i] = "blue4";
        }
        make_pdf(pos,edges,"view",color);
    };
    // DEBUG END
    print();
    auto a = 0;

    while(!found && !qs.empty() && !qt.empty()) {

        if(qs.size() < qt.size()) {
            // do one layer forward
            auto lastInLayer = qs.back();
            while(true) {
                auto v = qs.front(); qs.pop_front();
                space.insert(v);
                for(auto& e : adj[v]) {
                    edge_space.insert(minmax(v,e.to));
                    if(e.flow==e.cap || dist_s[e.to] <= dist_s[v]+1) continue;
                    dist_s[e.to] = dist_s[v]+1;
                    if(dist_t[e.to] < INF) found = true;
                    qs.push_back(e.to);
                }
                if(v==lastInLayer) break;
            }
        } else {
            // do one layer backward
            auto lastInLayer = qt.back();
            while(true) {
                auto v = qt.front(); qt.pop_front();
                space.insert(v);
                for(auto& e : adj[v]) {
                    edge_space.insert(minmax(v,e.to));
                    if(adj[e.to][e.rev].flow==adj[e.to][e.rev].cap || dist_t[e.to] <= dist_t[v]+1) continue;
                    dist_t[e.to] = dist_t[v]+1;
                    if(dist_s[e.to] < INF) found = true;
                    if(!found) qt.push_back(e.to);
                }
                if(v==lastInLayer) break;
            }
        }

        //print();
        auto a = 0;
    }

    // do rest forward
    while(found && !qs.empty()) {
        auto v = qs.front(); qs.pop_front();
        if(v==t||dist_t[v]==INF) continue;
        for(auto& e : adj[v]) {
            if(e.flow==e.cap || dist_t[e.to] >= dist_t[v]) continue;
            dist_s[e.to] = dist_s[v]+1;
            qs.push_back(e.to);
        }
    }

    cout << "bi\t" << space.size() << '\t' << edge_space.size() << endl;

    return dist_s;
}

int main() {
    ios::sync_with_stdio(false);

    int n = 500;
    int avg_deg = 5;
    int s = 33, t = 46;

    // girg params
    auto ple = 2.5;
    auto T = 0.2;
    int pseed = 0, wseed = 1, sseed = 2;
    pos = girgs::generatePositions(n,2,pseed,false);
    auto weights = girgs::generateWeights(n, ple, wseed, false);
    girgs::scaleWeights(weights,avg_deg,2,1/T);
    omp_set_num_threads(1);
    edges = girgs::generateEdges(weights,pos, 1/T, sseed);

    vector<vector<edge>> adj(n);
    for(auto [a,b] : edges) {
        adj[a].push_back(edge{a,b,0,1,adj[b].size()});
        adj[b].push_back(edge{b,a,0,1,adj[a].size()-1});
    }

    auto t1 = chrono::steady_clock::now();
    int flow =0;
    while(true) {
        auto dist = bfs_dist_bi2(s,t,n,adj);
        auto dist2= bfs_dist(s,t,n,adj);
        rep(i,n) assert(dist[i]==INF||dist2[i]==dist[i]);
        if(dist[t] == INF) break;
        vector<int> next(n,0);
        int aug;
        while(aug = dfs(s,INF,t,dist,next,adj)) flow += aug;
    }
    auto t2 = chrono::steady_clock::now();
    cout << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << endl;

    cout << flow << endl;

    return 0;
}
