
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
auto ghpairs(const string& file) { return baseFile(file) + ".ghpairs"; }

vector<pair<int,int>> readPairs(const string& pairsFile) {
    ifstream f(DATA_DIR + pairsFile);
    if(!f) cout << "error reading file " << pairsFile << endl, exit(0);
    int n; f>>n;
    vector<pair<int,int>> res(n);
    double dummy; // reads flow value
    for(auto& [s,t] : res) f>>s>>t>>dummy;
    return res;
}

vector<pair<int,int>> generatePairs(const vector<unsigned>& deg, int num, int mndeg, int mxdeg, mt19937& gen) {
    vector<int> candidates;
    rep(i,size(deg))
        if(mndeg<=deg[i] && deg[i]<=mxdeg)
            candidates.push_back(i);
    if(candidates.empty() || num > size(candidates)*(size(candidates)-1)/2) {
        cout << "not enough nodes in range " << mndeg << "-" << mxdeg << endl;
        cout << "candidates " << size(candidates) << endl;
        exit(0);
    }

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

// generate GH pairs
vector<pair<int,int>> generatePairs(const vector<pair<int,int>>& ghpairs, mt19937& gen, int num) {
    uniform_int_distribution dist(0ul, size(ghpairs)-1);
    set<size_t> choice;
    while(size(choice)<num) choice.insert(dist(gen));
    vector<pair<int,int>> res;
    for(auto idx : choice)
        res.emplace_back(ghpairs[idx].first, ghpairs[idx].second);
    return res;

}

// check if gh-pairs are present. if not compute them
void ensureGHpairs(const string& file) {
    {
        ifstream f(DATA_DIR + ghpairs(file));
        if(f) {
            cout << "existing gh-pairs found for " << file << endl;
            return;
        } else {
            cout << "no gh-pairs found for " << file << endl;
            cout << "computing new ones" << endl;
        }
    }

    auto graph = readGraph(file);

    vector<pair<int,int>> pairs;
    gomoryHuTree(graph, &pairs);

    Dinics::Dinics4Skip alg(graph);
    ofstream f(DATA_DIR + ghpairs(file));
    f << size(pairs) << '\n';
    for(auto [s,t] : pairs)
        f << s << '\t' << t << '\t' << alg.maxFlow(s,t) << '\n';
}

auto micros(chrono::steady_clock::time_point start, chrono::steady_clock::time_point end) {
    return chrono::duration_cast<chrono::microseconds>(end-start).count();
}

int main() {
    ios::sync_with_stdio(false);

    const vector<string> files {
            "fb-pages-tvshow.txt",    // 17K
            "girg10000.txt",          // 60K
            "soc-slashdot.txt",       // 360K
            "girg100000.txt",         // 600K
            "soc-flickr.txt",         // 3M
            "visualize-us.txt",       // 3M
            "dogster.txt",            // 8M
            "as-skitter.txt",         // 11M
            "actors.txt",             // 15M
            "brain.txt",              // 16M
    };

    const auto iters = 10; // number of measurements
    const auto seed = 1337;
    mt19937 gen(seed);

    ofstream log(DATA_DIR + "next_gen_exps.times"s);
    log << "inst,type,iter,alg,stage,time" << endl;
    for(const auto& file : files) {
        auto graph = readGraph(file);

        printGraphStats(graph,file);

        ensureGHpairs(file);
        auto allghpairs = readPairs(ghpairs(file));
        auto deg = degrees(graph);
        auto avg_deg = accumulate(all(deg), 0.0) / static_cast<double>(graph.n);
        auto [l_min, l_max] = pair(static_cast<int>(0.75*avg_deg), static_cast<int>(1.25*avg_deg));
        auto [h_min, h_max] = pair(static_cast<int>(10*avg_deg), static_cast<int>(100*avg_deg));
        cout << "low bounds     " << l_min << " / " << l_max << endl;
        cout << "high bounds    " << h_min << " / " << h_max << endl;

        array<vector<pair<int,int>>,3> pairs_by_type;
        array typeString{"low", "high", "gh"};
        pairs_by_type[0] = generatePairs(deg, iters, l_min, l_max, gen); // low
        pairs_by_type[1] = generatePairs(deg, iters, h_min, h_max, gen); // high
        pairs_by_type[2] = generatePairs(allghpairs, gen, iters); // gh

        // init all algs
        auto tt1 = chrono::steady_clock::now();
        Dinics::Dinics0Vanilla algUNI(graph);
        auto tt2 = chrono::steady_clock::now();
        Dinics::Dinics5OPT algBI(graph);
        auto tt3 = chrono::steady_clock::now();
        GoldbergTarjan::PushRelabel algPR(graph);
        auto tt4 = chrono::steady_clock::now();
        BoykovKolmogorov::BKAlgorithm algBK(graph);
        auto tt5 = chrono::steady_clock::now();

        //ofstream log(DATA_DIR + baseFile(file) + ".times");
        //log << "inst,type,iter,alg,stage,time" << endl;
        log << file << ",construct,0,Dinics,construct,"         << micros(tt1,tt2) << endl;
        log << file << ",construct,0,DinicsOPT,construct,"      << micros(tt2,tt3) << endl;
        log << file << ",construct,0,PushRelabel,construct,"    << micros(tt3,tt4) << endl;
        log << file << ",construct,0,BK-Algorithm,construct,"   << micros(tt4,tt5) << endl;

        array<vector<double>,3> answers;
        answers.fill(vector(iters,0.0));

        // DINICs
        rep(type,3) {
            ScopedTimer iterTimer("\tDinics "s + typeString[type]);
            rep(iter, iters) {
                auto [s,t] = pairs_by_type[type][iter];
                auto t1 = chrono::steady_clock::now();
                algUNI.resetFlow();
                auto t2 = chrono::steady_clock::now();
                answers[type][iter] = algUNI.maxFlow(s,t);
                auto t3 = chrono::steady_clock::now();
                log << file << ',' << typeString[type] << ',' << iter << ",Dinics,init," << micros(t1,t2) << '\n';
                log << file << ',' << typeString[type] << ',' << iter << ",Dinics,flow," << micros(t2,t3) << '\n';
            }
        }

        rep(type,3) {
            ScopedTimer iterTimer("\tDinicsOPT "s + typeString[type]);
            rep(iter, iters) {
                auto [s,t] = pairs_by_type[type][iter];
                auto t1 = chrono::steady_clock::now();
                algBI.resetFlow();
                auto t2 = chrono::steady_clock::now();
                auto f = algBI.maxFlow(s,t);
                if(f != answers[type][iter]) cout << "ERROR: wrong result\n", exit(0);
                auto t3 = chrono::steady_clock::now();
                log << file << ',' << typeString[type] << ',' << iter << ",DinicsOPT,init," << micros(t1,t2) << '\n';
                log << file << ',' << typeString[type] << ',' << iter << ",DinicsOPT,flow," << micros(t2,t3) << '\n';
            }
        }

        rep(type,3) {
            ScopedTimer iterTimer("\tPushRelabel "s + typeString[type]);
            rep(iter, iters) {
                auto [s,t] = pairs_by_type[type][iter];
                auto t1 = chrono::steady_clock::now();
                algPR.init(s,t);
                auto t2 = chrono::steady_clock::now();
                auto f= algPR.maxFlow(s,t);
                if(f != answers[type][iter]) cout << "ERROR: wrong result\n", exit(0);
                auto t3 = chrono::steady_clock::now();
                algPR.convertPreflow();
                auto t4 = chrono::steady_clock::now();
                algPR.deinit();
                log << file << ',' << typeString[type] << ',' << iter << ",PushRelabel,init," << micros(t1,t2) << '\n';
                log << file << ',' << typeString[type] << ',' << iter << ",PushRelabel,flow," << micros(t2,t3) << '\n';
                log << file << ',' << typeString[type] << ',' << iter << ",PushRelabel,convert," << micros(t3,t4) << '\n';
            }
        }

        rep(type,3) {
            ScopedTimer iterTimer("\tBK "s + typeString[type]);
            rep(iter, iters) {
                auto [s,t] = pairs_by_type[type][iter];
                auto t1 = chrono::steady_clock::now();
                auto f= algBK.maxFlow(s,t);
                if(f != answers[type][iter]) cout << "ERROR: wrong result\n", exit(0);
                auto t2 = chrono::steady_clock::now();
                log << file << ',' << typeString[type] << ',' << iter << ",BK-Algorithm,flow," << micros(t1,t2) << '\n';
            }
        }

    } // file

    return 0;
}
