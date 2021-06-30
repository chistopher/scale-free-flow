
#include <bits/stdc++.h>

#include <utils.h>
#include <Dinics.h>
#include <BoykovKolmogorov.h>
#include <PushRelabel.h>
#include <common.h>
#include <ScopedTimer.hpp>
#include <ghtree.h>
#include <ProgressIndicator.hpp>

#define rep(a, b)   for(int a = 0; a < (b); ++a)
#define all(a)      (a).begin(),(a).end()
//#define endl        '\n'

using namespace std;

// files that are passed are always relative to DATA_DIR
auto baseFile(const string& file) { return file.substr(0,size(file)-4); }
auto lpairs(const string& file) { return baseFile(file) + ".lpairs"; }
auto hpairs(const string& file) { return baseFile(file) + ".hpairs"; }
auto ghpairs(const string& file) { return baseFile(file) + ".ghpairs"; }
auto graph_file_from_pairs(const string& file) { return file.substr(0,file.rfind('.')) + ".txt"; }
// results are in .space

void writePairs(const vector<pair<int,int>>& pairs, const string& pairs_file) {
    auto graph = readGraph(graph_file_from_pairs(pairs_file));
    Dinics::Dinics4Skip alg(graph);
    ofstream f(DATA_DIR + pairs_file);
    f << size(pairs) << '\n';
    for(auto [s,t] : pairs)
        f << s << '\t' << t << '\t' << alg.maxFlow(s,t) << '\n';
}

vector<tuple<int,int,long long>> readPairs(const string& pairsFile) {
    ifstream f(DATA_DIR + pairsFile);
    if(!f) cout << "error reading file " << pairsFile << endl, exit(0);
    int n; f>>n;
    vector<tuple<int,int,long long>> res(n);
    for(auto& [s,t,w] : res) f>>s>>t>>w;
    return res;
}

vector<pair<int,int>> generatePairs(const vector<unsigned>& deg, int num, int mndeg, int mxdeg, mt19937& gen) {
    vector<int> candidates;
    rep(i,size(deg))
        if(mndeg<=deg[i] && deg[i]<=mxdeg)
            candidates.push_back(i);
    if(num > size(candidates)*(size(candidates)-1)/10)
        return {};

    uniform_int_distribution<int> dist(0, (int)size(candidates)-1);
    set<pair<int,int>> pairs;
    while(size(pairs)<num) {
        auto s = dist(gen);
        auto t = dist(gen);
        while (s == t) t = dist(gen);
        pairs.emplace(candidates[s],candidates[t]);
    }
    return {all(pairs)};
}

void writePairsLH(const string& file) {
    cout << "creating low and high deg pairs for " << file << endl;
    auto graph = readGraph(file);
    auto deg = degrees(graph);

    int seed = 1337;
    int num = 100;
    auto lowBounds = pair(1,10);
    auto highBounds = pair(100,5000);

    mt19937 gen(seed);
    auto lowDegPairs = generatePairs(deg, num, lowBounds.first, lowBounds.second, gen);
    auto highDegPairs = generatePairs(deg, num, highBounds.first, highBounds.second, gen);

    writePairs(lowDegPairs, lpairs(file));
    writePairs(highDegPairs, hpairs(file));
}

// get graph file as input
void writePairsGH(const string& file) {
    cout << "creating gh pairs file for " + file << endl;

    vector<pair<int,int>> pairs;
    gomoryHuTree(readGraph(file), &pairs);

    writePairs(pairs, ghpairs(file));
}


void computeSpace(const string& pairs_file) {
    cout << "measuring space for " << pairs_file << endl;
    auto graph = readGraph(graph_file_from_pairs(pairs_file));
    auto stw = readPairs(pairs_file);

    auto deg = degrees(graph);

    Dinics::Dinics0Vanilla algUNI(graph);
    Dinics::Dinics3Stamps algBI(graph);

    ofstream out(DATA_DIR + pairs_file + ".space.csv");
    out << "S,T,SDEG,TDEG,DIRECTION,ROUND," << Dinics::RoundStatistics::CSVHEADER << '\n';

    ProgressIndicator ind(size(stw));
    cout << "Compute space UNI" << endl;
    for (int i = 0; i < size(stw); ++i) {
        auto [s,t,w] = stw[i];
        auto f = algUNI.maxFlow(s,t,true);
        if(w!=f) cout << "mismatch " << s << ' ' << t << '\n', exit(0);
        int j = 0;
        for(auto& round : algUNI.stats())
            out << s << ',' << t << ',' << deg[s] << ',' << deg[t] << ",UNI," << j++ << ',' << round << '\n';
        algUNI.resetFlow();
        ind.tick(i);
    }

    cout << "Compute space BI" << endl;
    for (int i = 0; i < size(stw); ++i) {
        auto [s,t,w] = stw[i];
        auto f = algBI.maxFlow(s,t,true);
        if(w!=f) cout << "mismatch " << s << ' ' << t << '\n', exit(0);
        int j = 0;
        for(auto& round : algBI.stats())
            out << s << ',' << t << ',' << deg[s] << ',' << deg[t] << ",BI," << j++ << ',' << round << '\n';
        algBI.resetFlow();
        ind.tick(i);
    }
}

int main() {
    ios::sync_with_stdio(false);

    // lh pairs for a large graph
    const auto large = "as-skitter.txt";
    const auto medium = "soc-slashdot.txt";
    const auto girg = "girg10000.txt";

    //writePairsLH(large);
    //writePairsGH(medium);

    //computeSpace(lpairs(large));
    //computeSpace(hpairs(large));
    //computeSpace(ghpairs(medium));

    // timePairs(lpairs(large));

    auto graph_file = "as-skitter.txt";
    auto pairs_file = "as-skitter.pairs";
    auto graph = readGraph(graph_file);
    auto pairs = readPairs(pairs_file);
    Dinics::Dinics5OPT alg(graph);
    for(auto [s,t,w] : pairs) {
        if(alg.maxFlow(s,t) != w) {
            cout << "ERROR: pairs or alg are wrong" << endl;
            return 0;
        }
        alg.resetFlow();
    }

    computeSpace(pairs_file);

    return 0;
}
