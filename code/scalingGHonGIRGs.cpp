
#include <iostream>
#include <fstream>
#include <random>
#include <set>

#include <common.h>
#include <utils.h>
#include <ScopedTimer.hpp>
#include <Dinics.h>
#include <PushRelabel.h>
#include <BoykovKolmogorov.h>

#include <girgs/Generator.h>

#define rep(a, b)   for(int a = 0; a < (b); ++a)
#define all(a)      (a).begin(),(a).end()

using namespace std;

vector<pair<int,int>> generatePairs(const EdgeList& graph, int num, int mndeg, int mxdeg, mt19937& gen) {
    auto deg = degrees(graph);
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

auto micros(chrono::steady_clock::time_point start, chrono::steady_clock::time_point end) {
    return chrono::duration_cast<chrono::microseconds>(end-start).count();
}


int main() {
    ios::sync_with_stdio(false);

    int startSeed = 1337;
    int repetitions = 10; // number of girgs per n
    int num_pairs = 10; // number of s-t pairs per girg
    int avg_deg = 10;
    int min_n = 1'000;
    int max_n = 1'025'000;
    mt19937 gen(startSeed);

    auto sampleGirg = [avg_deg] (int n, int seed, int iter) {
        double ple = 2.8;
        double alpha = numeric_limits<double>::infinity(); // T=0
        int dimension = 1;
        int wseed = n + seed + iter + 12;
        int pseed = n + seed + iter + 123;
        int sseed = n + seed + iter + 1234; // irrelevant for T=0
        auto weights = girgs::generateWeights(n,ple,wseed,false);
        girgs::scaleWeights(weights, avg_deg, dimension, alpha);
        const auto positions = girgs::generatePositions(n,dimension,pseed,false);
        EdgeList res;
        res.edges = girgs::generateEdges(weights,positions,alpha,sseed);
        res.n = n;
        return res;
    };

    ofstream f(DATA_DIR + string("other/scaling.times"));
    f << "n,iter,alg,time\n";
    for (int n = min_n; n < max_n; n*=2) {
        cout << "n " << n << endl;
        for(int iter=0; iter<repetitions; ++iter) {
            ScopedTimer timer("\titer " + to_string(iter));

            auto graph = sampleGirg(n, startSeed, iter);
            Dinics::Dinics0Vanilla  din0(graph);
            Dinics::Dinics1Bi       din1(graph);
            Dinics::Dinics2Reset    din2(graph);
            Dinics::Dinics3Stamps   din3(graph);
            Dinics::Dinics4Skip     din4(graph);
            Dinics::Dinics5OPT      din5(graph);
            GoldbergTarjan::PushRelabel algPR(graph);
            BoykovKolmogorov::BKAlgorithm algBK(graph);

            auto pairs = generatePairs(graph, num_pairs, 10, 20, gen);

            // time all dinics
            auto timeAlg = [&](auto& alg, const string& name) {
                auto elp = 0ll;
                for(auto [s,t] : pairs) {
                    alg.resetFlow();
                    auto t1 = chrono::steady_clock::now();
                    alg.maxFlow(s, t);
                    auto t2 = chrono::steady_clock::now();
                    elp += micros(t1,t2);
                }
                f << n << ',' << iter << ',' << name << ',' <<  elp / num_pairs << endl;
            };

            timeAlg(din0, "Dinics0Vanilla");
            timeAlg(din1, "Dinics1Bi");
            timeAlg(din2, "Dinics2Reset");
            timeAlg(din3, "Dinics3Stamp");
            timeAlg(din4, "Dinics4Skip");
            timeAlg(din5, "Dinics5OPT");

            // time PR
            {
                auto elp = 0ll;
                for(auto [s,t] : pairs) {
                    algPR.init(s,t);
                    auto t1 = chrono::steady_clock::now();
                    algPR.maxFlow(s, t);
                    auto t2 = chrono::steady_clock::now();
                    elp += micros(t1,t2);
                    algPR.deinit();
                }
                f << n << ',' << iter << ',' << "PushRelabel" << ',' <<  elp / num_pairs << endl;
            }

            // time BK
            {
                auto elp = 0ll;
                for(auto [s,t] : pairs) {
                    auto t1 = chrono::steady_clock::now();
                    algBK.maxFlow(s, t);
                    auto t2 = chrono::steady_clock::now();
                    elp += micros(t1,t2);
                }
                f << n << ',' << iter << ',' << "BK-Algorithm" << ',' <<  elp / num_pairs << endl;
            }
        }
    }


    return 0;
}
