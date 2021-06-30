
#include <bits/stdc++.h>

#include <utils.h>
#include <Dinics.h>
#include <common.h>
#include <ScopedTimer.hpp>
#include <ProgressIndicator.hpp>

#define rep(a, b)   for(int a = 0; a < (b); ++a)
#define all(a)      (a).begin(),(a).end()
//#define endl        '\n'

using namespace std;

void writePairs(const EdgeList& graph, const vector<pair<int,int>>& pairs, const string& pairs_file) {
    Dinics::Dinics4Skip alg(graph);
    ofstream f(DATA_DIR + pairs_file);
    f << size(pairs) << '\n';
    for(auto [s,t] : pairs) {
        f << s << '\t' << t << '\t' << alg.maxFlow(s,t) << '\n';
        alg.resetFlow();
    }
}

// this is for the benchmarks in the table
// we also write the pairs to file to compute space stats with the 'pairs' executable later
// result file: data/as-skitter.pairs
int main() {
    ios::sync_with_stdio(false);

    const auto file = "as-skitter.txt";
    auto graph = readGraph(file);
    //printGraphStats(graph);

    auto deg_mn = 10;
    auto deg_mx = 16;
    auto seed = 2500;

    mt19937 gen(seed);

    auto pairs = generatePairs(graph, 1000, deg_mn, deg_mx, gen);
    //writePairs(graph, pairs, "as-skitter.pairs");

    Dinics::Dinics0Vanilla alg(graph);
    //Dinics::Dinics1Bi alg(graph);
    //Dinics::Dinics2Reset alg(graph);
    //Dinics::Dinics3Stamps alg(graph);
    //Dinics::Dinics4Skip alg(graph);
    //Dinics::Dinics5OPT alg(graph);

    for(auto [s,t] : pairs) {
        alg.resetFlow();
        alg.maxFlow(s,t);
    }

    return 0;
}
