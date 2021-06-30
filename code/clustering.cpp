
#include <bits/stdc++.h>

#include <utils.h>
#include <ScopedTimer.hpp>
#include <ProgressIndicator.hpp>
#include <common.h>
#include <Dinics.h>
#include <PushRelabel.h>

#include <girgs/Generator.h>

#define rep(a, b)   for(int a = 0; a < (b); ++a)
#define all(a)      (a).begin(),(a).end()
//#define endl        '\n'

using namespace std;

// each edge (u,v) costs 1/deg(u) + 1/deg(v)
// TODO handle directed edges differently?
void normalizedWeights(EdgeList& graph) {
    auto deg = degrees(graph);
    if(!graph.hasWeights) {
        graph.hasWeights = true;
        graph.weights.resize(size(graph.edges));
    }
    rep(i,size(graph.edges)) {
        auto [u,v] = graph.edges[i];
        graph.weights[i] = 1.0 / deg[u] + 1.0 / deg[v];
    }
}

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
            if(!reverse && e.isSaturated()) continue;
            if(reverse && e.rev->isSaturated()) continue;
            visited[e.to] = true;
            reach.push_back(e.to);
        }
    }
    return reach;
}

// overload for Push Relabel Network
template<>
vector<int> residualReachable<GoldbergTarjan::LinearizedNetwork>(const GoldbergTarjan::LinearizedNetwork& g, int s_or_t, bool reverse) {
    vector visited(g.n, false);
    visited[s_or_t] = true;
    int qB = 0; // queue begin
    vector<int> reach{s_or_t};
    while(qB<reach.size()) {
        auto v = reach[qB++];
        for(auto [eb, ee] = g.adj[v]; eb!=ee; ++eb) {
            auto& e = g.data[eb];
            if(visited[e.to]) continue;
            if(!reverse && e.residual_capacity<1e-6) continue;
            if(reverse && g.data[e.reverse_edge].residual_capacity<1e-6) continue;
            visited[e.to] = true;
            reach.push_back(e.to);
        }
    }
    return reach;
}

// bigger values means more clusters
vector<int> clustering(const EdgeList& original, double alpha, bool writeFiles = false) {
    // output stuff
    vector<pair<int,int>> pairs;

    // add super sink
    int superSink = static_cast<int>(original.n);
    auto graph = original;
    graph.n++;
    if(!original.hasWeights) {
        graph.hasWeights = true;
        graph.weights.resize(size(graph.edges), 1.0);
    }
    rep(v,original.n) {
        graph.edges.emplace_back(v,superSink);
        graph.weights.push_back(alpha);
    }
    // revert edges and weights so that edges to t are first in adj list
    reverse(all(graph.edges));
    reverse(all(graph.weights));

    // sort nodes desc by degree to improve performance
    vector<int> nodes(original.n);
    iota(all(nodes),0);
    auto deg = degrees(graph);
    sort(all(nodes),[&](int a, int b) { return deg[a] > deg[b]; });

    // GUSFILED GH TREE alg trimmed down
    int skipped = 0;
    vector<int> parent(original.n, superSink); // parent in GH tree
    Dinics::Dinics4Skip alg(graph);
    //GoldbergTarjan::PushRelabel alg(graph);
    ProgressIndicator ind(original.n); int done = 0;
    for(auto s : nodes) {
        ind.tick(done++);
        if(parent[s] != superSink) {
            skipped++;
            continue;
        }

        alg.resetFlow();
        alg.maxFlow(s,superSink);
        pairs.emplace_back(s,superSink);

        for(auto v : residualReachable(alg.graph(), s)) {
            if(parent[v]==superSink && v!=s)
                parent[v] = s;
        }
    }
    cout << "skipped " << skipped << endl;

    // write stuff
    if(writeFiles) {
        exportDIMACS(graph, "other/clustering.max");
        ofstream pairsFile(DATA_DIR + "other/clustering.pairs"s);
        pairsFile << size(pairs) << '\n';
        for(auto [a,b] : pairs) pairsFile << a << ' ' << b << '\n';
    }

    // hang all nodes directly below their representative
    rep(i,original.n) if(parent[i] == superSink) parent[i] = i;
    auto findRep = [&] (auto&& f, int i) {
        if(parent[i]==i) return i;
        return parent[i] = f(f,parent[i]);
    };
    rep(i,original.n) findRep(findRep,i);

    return parent;
}

// returns the next clustering with at least number clusters
vector<int> binarySearchClusters(const EdgeList& graph, int number) {
    ScopedTimer timer("compute clustering");
    double low = 0.0; // too few clusters
    double high = 1.0; // too many clusters
    while(high-low > 1e-5) {
        auto mid = (high+low)/2;
        auto cls = clustering(graph, mid);
        sort(all(cls));
        cls.erase(unique(all(cls)),end(cls));
        cout << low << '\t' << high << '\t' << size(cls) << endl;
        if(size(cls) < number)
            low = mid;
        else
            high = mid;
    }
    return clustering(graph, high);
}

auto generateGIRG() {
    int n = 5'000;
    int deg = 5;
    int d = 2;
    double ple = 2.9;
    double alpha = numeric_limits<double>::infinity();
    int pseed = 1336;
    int wseed = 5;
    int eseed = 1339;
    auto positions = girgs::generatePositions(n,d,pseed,false);
    auto weights = girgs::generateWeights(n,ple,wseed,false);
    girgs::scaleWeights(weights, deg, d, alpha);
    auto edges = girgs::generateEdges(weights,positions,alpha,eseed);

    EdgeList g;
    g.n = n;
    g.edges = edges;
    return pair(g, positions);
}

void printClusterStats(const EdgeList& graph, const vector<int>& cluster) {
    auto deg = degrees(graph);
    map<int,int> clusterSize;
    for(auto r : cluster) clusterSize[r]++;
    auto zeros = count(all(deg), 0);
    cout << "clusters           " << size(clusterSize) << endl;
    cout << "unconnected nodes  " << zeros << endl;
    cout << "largest clusters   ";
    vector<int> cls;
    for(auto kv : clusterSize) cls.push_back(kv.second);
    sort(all(cls),greater<>());
    rep(i,min<int>(10,(int)size(cls))) cout << cls[i] << ' ';
    cout << endl;
    cout << "avg cluster size   " << accumulate(all(cls),0.0)/size(cls) << endl;
    cout << "total cluster size " << accumulate(all(cls),0ll) << endl;
    for(int i = 1; i<8; ++i) cout << "size " << i << " clusters    " << count(all(cls), i) << endl;
}

int main() {
    ios::sync_with_stdio(false);

    const string file = "other/clustering.max";

    auto [graph, pos] = generateGIRG();
    map<int,int> ids;
    rep(i,graph.n) ids[i] = i;

    {
        // print stats
        cout << "===========================" << endl;
        cout << "n       " << graph.n << endl;
        cout << "m       " << size(graph.edges) << endl;
        cout << "avg deg " << 2.0 * static_cast<double>(size(graph.edges)) / graph.n << endl;
        cout << "===========================" << endl;
        graph = largestCC(graph, &ids);
        cout << "n       " << graph.n << endl;
        cout << "m       " << size(graph.edges) << endl;
        cout << "avg deg " << 2.0 * static_cast<double>(size(graph.edges)) / graph.n << endl;
        cout << "===========================" << endl;
    }

    map<int,int> ids_inv;
    for(auto [k,v] : ids) ids_inv[v] = k;

    //normalizedWeights(graph);

    //auto cluster = binarySearchClusters(graph, 12);
    auto cluster = clustering(graph, 0.02,true);

    bool debug = false;
    if(debug) {
        printClusterStats(graph, cluster);

        map<int,int> colors, sz;
        for(auto v : cluster) sz[v]++;
        vector<pair<int,int>> cls(begin(sz), end(sz));
        sort(all(cls),[](auto p1, auto p2) { return p1.second > p2.second; });
        for(auto [c,csz] : cls) {
            if(colors.find(c)==colors.end())
                colors[c] = colors.size() < 12 ? colors.size()+1 : 0;
        }
        ofstream f("graph.dot");
        f<<"graph {\nscale=3000;\n";
        rep(v,graph.n) f << v << " [ "
                         << "cluster=" << cluster[v]
                         << " style =\"filled\" fillcolor=\"/paired12/" << colors[cluster[v]] << '\"'
                         << " pos=\"" << pos[ids_inv[v]][0] << ',' << pos[ids_inv[v]][1] << "\"]\n";
        for(auto [a,b] : graph.edges) f << a << " -- " << b << '\n';
        f << '}';
        f.close();

        string command = "neato -n -Tpdf -Nlabel='' graph.dot -o graph.pdf";
        system(command.c_str());
    }

    return 0;
}
