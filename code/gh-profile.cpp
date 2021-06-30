
#include <bits/stdc++.h>
#include <cstdlib>

#include <ProgressIndicator.hpp>
#include <utils.h>
#include <Dinics.h>
#include <PushRelabel.h>
#include <ScopedTimer.hpp>


#define rep(a, b)   for(int a = 0; a < (b); ++a)
#define all(a)      (a).begin(),(a).end()
//#define endl        '\n'

using namespace std;

// impl for Dinics
// TODO this does only work for directed because of the adj-rev opt
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

// impl for PR
template<>
vector<int> residualReachable(const GoldbergTarjan::LinearizedNetwork& g, int s_or_t, bool reverse) {
    vector visited(g.n, false);
    visited[s_or_t] = true;
    int qB = 0; // queue begin
    vector<int> reach{s_or_t};
    while(qB<reach.size()) {
        auto v = reach[qB++];
        auto [eb, ee] = g.adj[v];
        for(auto i = eb; i<ee; ++i) {
            GoldbergTarjan::Edge e = g.data[i];
            if(visited[e.to]) continue;
            if(!reverse && e.residual_capacity<=0) continue;
            auto rev_e = g.data[e.reverse_edge];
            if(reverse && rev_e.residual_capacity<=0) continue;
            visited[e.to] = true;
            reach.push_back(e.to);
        }
    }
    return reach;
}

template<class CutOracle>
auto gomoryHuTreeMeasuringProxy(CutOracle& f, int n) {
    vector parent(n, 0);
    for (int s = 1; s < n; ++s) {
        auto t = parent[s];
        auto s_side_cut = f(s,t);
        for(auto v : s_side_cut) { // go through S side of cut
            // hang all siblings of s below s
            if(v!=s && parent[v]==t)
                parent[v] = s;
            // hang t below s and s below ts former parent
            if(v==parent[t]) {
                parent[s] = parent[t];
                parent[t] = s;
            }
        }
    }
    return parent;
}

// returns all 0 <= i < n that are not in vec
vector<int> complement(const vector<int>& vec, int n) {
    assert(all_of(all(vec), [n](int v){ return 0<v&&v<n; }));
    vector in_vec(n, false);
    for(auto v : vec) in_vec[v] = true;
    vector<int> res;
    rep(i,n) if(!in_vec[i]) res.push_back(i);
    return res;
}

auto cutValue(const EdgeList& g, const vector<int>& s_side) {
    vector in_s(g.n, false);
    for(auto v : s_side) in_s[v] = true;
    auto res = 0.0;
    rep(i,size(g.edges)) {
        auto w = g.hasWeights ? g.weights[i] : 1.0;
        auto [u,v] = g.edges[i];
        if((g.isDirected && in_s[u] && !in_s[v]) ||
           (!g.isDirected && in_s[u] != in_s[v]))
            res += w;
    }
    return res;
}

auto micros(chrono::steady_clock::time_point start, chrono::steady_clock::time_point end) {
    return chrono::duration_cast<chrono::microseconds>(end-start).count();
}

// used to benchmark different cut-id versions for PR
// results are on stdout
int main() {
    ios::sync_with_stdio(false);

    const string file = "soc-slashdot.txt";
    //const string file = "as-skitter.txt";
    //const string file = "girg100000.txt";
    //const string file = "girglowdeg.txt";
    //const string file = "fb-pages-tvshow.txt";

    auto graph = readGraph(file);

    {
        // print stats
        cout << "===========================" << endl;
        cout << "file    " << file << endl;
        cout << "n       " << graph.n << endl;
        cout << "m       " << size(graph.edges) << endl;
        cout << "avg deg " << 2.0 * static_cast<double>(size(graph.edges)) / graph.n << endl;
        cout << "===========================" << endl;
    }

    const auto n = graph.n;

    // Dinics validator functions
    Dinics::Dinics5OPT alg(graph);
    auto s_side_from_s = [&alg](int s, int t){
        alg.resetFlow();
        alg.maxFlow(s,t);
        return residualReachable(alg.graph(), s);
    };
    auto t_side_from_s = [&alg, n](int s, int t){
        alg.resetFlow();
        alg.maxFlow(s,t);
        return complement(residualReachable(alg.graph(), t, true), n);
    };
    auto s_side_from_t = [&alg](int s, int t){
        alg.resetFlow();
        alg.maxFlow(t,s);
        return residualReachable(alg.graph(), s, true);
    };
    auto t_side_from_t = [&alg, n](int s, int t){
        alg.resetFlow();
        alg.maxFlow(t,s);
        return complement(residualReachable(alg.graph(), t), n);
    };

    vector<int> parents, parents_rev; // gh trees when taking sside/tside cut each time
    {
        ScopedTimer timer("DINICS best");
        parents = gomoryHuTreeMeasuringProxy(s_side_from_s, n);
    }
    {
        ScopedTimer timer("DINICS rev");
        parents_rev = gomoryHuTreeMeasuringProxy(t_side_from_s, n);
    }


    // time accumulators
    auto init = 0ll;
    auto flow = 0ll;
    auto convert = 0ll;
    auto cut = 0ll;

    auto printStats = [&] {
        cout << "init    " << init/1000 << endl;
        cout << "flow    " << flow/1000 << endl;
        cout << "convert " << convert/1000 << endl;
        cout << "cut     " << cut/1000 << endl;
        init = flow = convert = cut = 0;
    };

    // same for PR
    GoldbergTarjan::PushRelabel algPR(graph);
    // CONVERT
    auto PR_s_side_from_s = [&](int s, int t){
        auto tt1 = chrono::steady_clock::now();
        algPR.init(s,t);
        auto tt2 = chrono::steady_clock::now();
        algPR.maxFlow(s,t);
        auto tt3 = chrono::steady_clock::now();
        algPR.convertPreflow();
        auto tt4 = chrono::steady_clock::now();
        auto res = residualReachable(algPR.graph(), s);
        auto tt5 = chrono::steady_clock::now();
        algPR.deinit();

        init += micros(tt1, tt2);
        flow += micros(tt2, tt3);
        convert += micros(tt3, tt4);
        cut += micros(tt4, tt5);

        return res;
    };

    {
        ScopedTimer timer("PR CONVERT");
        auto res = gomoryHuTreeMeasuringProxy(PR_s_side_from_s, n);
        if(res!=parents) cout << "wrong tree" << endl, exit(0);
    }
    printStats();

    // T-SIDE
    auto PR_t_side_from_s = [&](int s, int t){
        auto tt1 = chrono::steady_clock::now();
        algPR.init(s,t);
        auto tt2 = chrono::steady_clock::now();
        algPR.maxFlow(s,t);
        auto tt3 = chrono::steady_clock::now();
        // algPR.convertPreflow();
        auto tt4 = chrono::steady_clock::now();
        auto res = complement(residualReachable(algPR.graph(), t, true), n);
        auto tt5 = chrono::steady_clock::now();
        algPR.deinit();

        init += micros(tt1, tt2);
        flow += micros(tt2, tt3);
        // convert += micros(tt3, tt4);
        cut += micros(tt4, tt5);

        return res;
    };

    {
        ScopedTimer timer("PR T-SIDE");
        auto res = gomoryHuTreeMeasuringProxy(PR_t_side_from_s, n);
        if(res!=parents_rev) cout << "wrong tree" << endl, exit(0);
    }
    printStats();


    // SWAP
    auto PR_s_side_from_t = [&](int s, int t){
        auto tt1 = chrono::steady_clock::now();
        algPR.init(t,s);
        auto tt2 = chrono::steady_clock::now();
        algPR.maxFlow(t,s);
        auto tt3 = chrono::steady_clock::now();
        // algPR.convertPreflow();
        auto tt4 = chrono::steady_clock::now();
        auto res = residualReachable(algPR.graph(), s, true);
        auto tt5 = chrono::steady_clock::now();
        algPR.deinit();

        init += micros(tt1, tt2);
        flow += micros(tt2, tt3);
        // convert += micros(tt3, tt4);
        cut += micros(tt4, tt5);

        //algPR.init(t,s);
        //algPR.maxFlow(t,s);
        //algPR.deinit();
        //return residualReachable(algPR.graph(), s, true);
        return res;
    };
    {
        ScopedTimer timer("PR SWAP");
        auto res = gomoryHuTreeMeasuringProxy(PR_s_side_from_t, n);
        if(res!=parents) cout << "wrong tree" << endl, exit(0);
    }
    printStats();

    // this is useless
    auto PR_t_side_from_t = [&algPR, n](int s, int t){
        // needs swap and conversion
        algPR.init(t,s);
        algPR.maxFlow(t,s);
        algPR.convertPreflow();
        algPR.deinit();
        return complement(residualReachable(algPR.graph(), t), n);
    };


    // testing
    /*
    mt19937 gen(1234235);
    uniform_int_distribution dist(0,(int)n-1);
    rep(i,100000) {
        auto s = dist(gen);
        auto t = dist(gen);
        while(s==t) t = dist(gen);

        alg.resetFlow();
        auto f = alg.maxFlow(s,t);

        auto p1 = s_side_from_s(s,t);
        auto p2 = t_side_from_s(s,t);
        auto p3 = s_side_from_t(s,t);
        auto p4 = t_side_from_t(s,t);
        sort(all(p1)); sort(all(p2));
        sort(all(p3)); sort(all(p4));

        cout << "=== " << s << " " << t << " ===" << endl;
        cout << (p1 == p2) << endl;
        cout << (p1 == p3) << endl;
        cout << (p2 == p4) << endl;

        cout << "vals: " << endl;
        cout << f << endl;
        for(const auto& cut : {p1,p2,p3,p4})
            cout << cutValue(graph, cut) << " (" << size(cut) << ")\n";

        auto q1 = PR_s_side_from_s(s,t);
        auto q2 = PR_t_side_from_s(s,t);
        auto q3 = PR_s_side_from_t(s,t);
        auto q4 = PR_t_side_from_t(s,t);
        sort(all(q1)); sort(all(q2));
        sort(all(q3)); sort(all(q4));

        cout << "PR check ... " << endl;
        vector<vector<int>> parr{p1,p2,p3,p4};
        vector<vector<int>> qarr{q1,q2,q3,q4};
        rep(i,4) {
            auto& p = parr[i];
            auto& q = qarr[i];
            if(p==q) continue;
            cout << "PR differs" << endl;
            for(auto cut : {p,q})
                cout << '\t' << cutValue(graph, cut) << " (" << size(cut) << ")\n";
        }
    }
    */


    return 0;
}
