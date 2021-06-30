
#include <bits/stdc++.h>

#include <ScopedTimer.hpp>
#include <ProgressIndicator.hpp>
#include <common.h>
#include <utils.h>
#include <Dinics.h>
#include <BoykovKolmogorov.h>
#include <PushRelabel.h>

#define rep(a, b)   for(int a = 0; a < (b); ++a)
#define all(a)      (a).begin(),(a).end()
//#define endl        '\n'

using namespace std;

auto baseFile(const string& file) { return file.substr(0,size(file)-4); }

EdgeList ErdosRenyi(size_t n, double p, int seed = 13) {
    assert(0<p && p<1.0);
    EdgeList g;
    g.n = n;
    g.isDirected = true;
    mt19937 gen(seed);
    geometric_distribution<size_t> dist(p);
    auto e = dist(gen);
    while(e<n*n) {
        g.edges.emplace_back(e/n, e%n);
        e += 1 + dist(gen);
    }
    return g;
}

vector<pair<int,int>> randomPairs(int n, int num, mt19937& gen) {
    set<pair<int,int>> res;
    uniform_int_distribution dist(0, n-1);
    while(res.size()<num) {
        auto a = dist(gen);
        auto b = dist(gen);
        while(a==b) b = dist(gen);
        res.emplace(a,b);
    }
    return {all(res)};
}

vector<pair<int,int>> pairsFromFile(const string& pairsFile) {
    ifstream f(DATA_DIR + pairsFile);
    if(!f) cout << "cannot read file " << pairsFile << endl, exit(0);
    int num; f>>num;
    vector<pair<int,int>> res(num);
    for(auto& [s,t] : res) f>>s>>t;
    return res;
}

// generates a network from Ahuja1997
auto syntheticNetworkLayered(int W, int L, int d, int seed = 14) {
    mt19937 gen(seed);

    EdgeList g;
    g.n = W*L+2;
    g.hasWeights = true;
    g.isDirected = true;

    auto source = g.n-2;
    auto sink = g.n-1;
    auto INF = 1e9;
    uniform_int_distribution capDist(500,10'000);
    uniform_int_distribution degDist(1,2*d-1);
    uniform_int_distribution neiDist(0,W-1);
    auto edge = [&](int a, int b, bool infCap = false) {
        g.edges.emplace_back(a,b);
        g.weights.push_back(infCap ? INF : capDist(gen));
    };

    rep(l,L) {
        rep(v,W) {
            if(l==L-1) { edge(l*W+v, sink, true); continue; }
            if(l==0) edge(source, l*W+v, true);
            auto deg = degDist(gen);
            rep(i,deg) edge(l*W+v, (l+1)*W+neiDist(gen));
        }
    }

    return g;
}

auto addRandomWeights(const EdgeList& g, int seed = 15) {
    mt19937 gen(seed);
    uniform_int_distribution dist(500,10'000);
    auto g2 = g;
    g2.hasWeights = true;
    g2.weights.resize(g2.edges.size());
    for(auto& w : g2.weights) w = dist(gen);
    return g2;
}

// adds two more nodes that are connection to size many nodes in the given graph
// if graph is directed, edges will only be from n, and to n+1
auto addSuperTerminals(const EdgeList& g, int size, int seed = 16) {
    mt19937 gen(seed);
    uniform_int_distribution dist(0, (int)g.n-1);
    set<int> superS, superT;
    while(superS.size() < size) superS.insert(dist(gen));
    while(superT.size() < size) {
        auto cand = dist(gen);
        if(superS.find(cand)==superS.end()) superT.insert(cand);
    }
    auto g2 = g;
    auto s = g.n;
    auto t = g.n+1;
    g2.n += 2;
    auto INF = 1e9;
    for(auto v : superS) {
        g2.edges.emplace_back(s,v);
        if(g2.hasWeights) g2.weights.push_back(INF);
    }
    for(auto v : superT) {
        g2.edges.emplace_back(v,t);
        if(g2.hasWeights) g2.weights.push_back(INF);
    }
    return g2;
}

void generateNetworks() {
    // same params as in Ahuja1997
    exportDIMACS(syntheticNetworkLayered(71,141,10), "other/layered10000.max", 71*141, 71*141+1);
    // some variations on ER
    auto ER_SIZE = 100'000;
    auto ER_P = 0.002;
    auto ER_NAME = "other/er" + to_string(ER_SIZE);
    auto ER_GRAPH = ErdosRenyi(ER_SIZE, ER_P);
    exportDIMACS(ER_GRAPH, ER_NAME + ".max");
    exportDIMACS(addRandomWeights(ER_GRAPH), ER_NAME + "_weighted.max");
    exportDIMACS(addSuperTerminals(ER_GRAPH, 2000), ER_NAME + "_super.max", ER_SIZE, ER_SIZE+1); // s,t are last two vertices
}

auto micros(chrono::steady_clock::time_point start, chrono::steady_clock::time_point end) {
    return chrono::duration_cast<chrono::microseconds>(end-start).count();
}

// run a bunch of terminal pairs
template<class ALG>
void computeMultiple(ALG& alg, const vector<pair<int,int>>& terminals) {
    for(auto [s,t] : terminals) {
        alg.resetFlow();
        alg.maxFlow(s,t);
    }
}

void computeMultiple(BoykovKolmogorov::BKAlgorithm& alg, const vector<pair<int,int>>& terminals) {
    for(auto [s,t] : terminals) alg.maxFlow(s,t);
}

void computeMultiple(GoldbergTarjan::PushRelabel& alg, const vector<pair<int,int>>& terminals) {
    for(auto [s,t] : terminals) {
        alg.init(s,t);
        alg.maxFlow(s,t);
        alg.convertPreflow();
        alg.deinit();
    }
}

int main() {
    ios::sync_with_stdio(false);

    //generateNetworks();

    vector<string> files {
            "other/er100000.max",
            "other/er100000_weighted.max",
            "other/er100000_super.max",
            "other/layered10000.max",
            "other/roadNet-PA.max",
            "other/liver.n6c100.max", // solution: 625447 (tested)
            "other/clustering.max",
    };

    const auto iters = 50; // number of s-t measurements
    auto clusteringPairs = pairsFromFile("other/clustering.pairs");
    const auto seed = 1337;
    mt19937 gen(seed);

    ofstream log(DATA_DIR + "other/combined.times2"s);
    log << "inst,iter,alg,time" << endl;
    for(auto file : files) {
        auto [graph,s,t] = importDIMACS(file);
        printGraphStats(graph, file);
        auto isClustering = (file == "other/clustering.max");

        // pairs for each iteration can be given or random
        auto pairs = s!=-1 ? vector(iters, pair(s,t)) : randomPairs(graph.n, iters, gen);

        for (int iter = 0; iter < iters; ++iter) {
            cout << "\tstarting iteration " << iter << endl;
            ScopedTimer iterTimer("\titeration " + to_string(iter));

            auto [source,sink] = pairs[iter];

            Dinics::Dinics4Skip algOPT(graph);
            Dinics::Dinics0Vanilla algUNI(graph);
            GoldbergTarjan::PushRelabel algPR(graph);
            if(!isClustering) algPR.init(source, sink);
            BoykovKolmogorov::BKAlgorithm algBK(graph);

            cout << "\t\tinit all algs done in " << iterTimer.elapsed() << endl;

            auto f = algOPT.maxFlow(source,sink); algOPT.resetFlow();
            auto timeAlg = [&](auto &alg, const string& name) {
                ScopedTimer timer("\t\ttiming alg " + name);
                auto t1 = chrono::steady_clock::now();
                double res = 0;
                if(!isClustering) res = alg.maxFlow(source, sink);
                else computeMultiple(alg, clusteringPairs);
                auto t2 = chrono::steady_clock::now();
                auto eps = micros(t1,t2);
                log << baseFile(file).substr(6) << ',' << iter << ',' << name << ',' << eps << endl;
                if(!isClustering && f!=res) {
                    cout << "WARNING: flow values differ for alg " << name << endl;
                    cout << "\tResult: " << f << endl;
                    cout << "\tGiven:  " << res << endl;
                }
            };

            timeAlg(algUNI, "Dinics");
            timeAlg(algOPT, "DinicsOPT");
            timeAlg(algPR, "PushRelabel");
            timeAlg(algBK, "BK-Algorithm");
        }
    }

    return 0;
}
