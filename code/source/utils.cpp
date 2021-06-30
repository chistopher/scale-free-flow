
#include <utils.h>

#include <fstream>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <set>

#include <common.h>
#include <ScopedTimer.hpp>

using namespace std;

// header is n m hasWeights(0 or 1) isDirected(0 or 1)
EdgeList readGraph(const string& file) {
    ifstream f(DATA_DIR + file);
    if(!f) cout << "\nerror reading file " << file << endl, exit(0);
    char* data;
    {
        f.seekg(0, ios::end);;
        size_t length = f.tellg();
        f.seekg(0, ios::beg);
        data = (char*)malloc(length);
        f.read(data,length);
    }

    auto readInt = [c=data](auto &number) mutable {
        while(*c == ' ' || *c == '\n' || *c == '\t') ++c;
        number = 0;
        for (; '0'<=*c && *c<='9'; ++c)
            number = number*10 + *c - '0';
    };

    int n, m, hasW, isDir;
    readInt(n);
    readInt(m);
    readInt(hasW);
    readInt(isDir);

    EdgeList graph;
    graph.n = n;
    graph.isDirected = isDir;
    graph.hasWeights = hasW;
    graph.edges.resize(m);
    graph.weights.resize(hasW ? m : 0);

    for (int i = 0; i < m; ++i) {
        readInt(graph.edges[i].first);
        readInt(graph.edges[i].second);
        if(hasW) readInt(graph.weights[i]);
    }
    free(data);
    return graph;
}

std::vector<unsigned> degrees(const EdgeList& edges) {
    vector<unsigned> deg(edges.n, 0);
    for(auto [a,b] : edges.edges) {
        assert(0<=a && a<edges.n);
        assert(0<=b && b<edges.n);
        deg[a]++;
        deg[b] += !edges.isDirected;
    }
    return deg;
}

EdgeList largestCC(const EdgeList &graph, std::map<int,int>* out) {
    // TODO use linearized network
    auto n = graph.n;
    vector<vector<int>> g(n);
    for(auto [a,b] : graph.edges) {
        g[a].push_back(b);
        if(!graph.isDirected) g[b].push_back(a);
    }

    // find all connected components
    auto maxComp = pair(1,0); // (size,rep)
    vector<int> comp(n,-1);
    for(int v=0; v<n; ++v) {
        if(comp[v]!=-1) continue;
        comp[v] = v;
        int compSize = 1;
        vector<int> stack{v};
        while(!stack.empty()) {
            auto u = stack.back(); stack.pop_back();
            for(auto nei : g[u]) {
                assert(comp[nei]==-1 || comp[nei]==v);
                if(comp[nei]==v) continue;
                comp[nei] = v;
                stack.push_back(nei);
                compSize++;
            }
        }
        maxComp = max(maxComp, pair(compSize,v));
    }

    // filter edges and assign new ids
    EdgeList res;
    res.n = static_cast<size_t>(maxComp.first);
    res.hasWeights = graph.hasWeights;
    res.isDirected = graph.isDirected;

    map<int,int> ids;
    for(int i=0; i<size(graph.edges); ++i) {
        auto [a,b] = graph.edges[i];
        assert(comp[a]!=-1);
        assert(comp[a]==comp[b]);
        if(comp[a] != maxComp.second) continue;
        for(int v : {a,b})
            if(ids.find(v) == ids.end())
                ids[v] = static_cast<int>(size(ids));
        res.edges.emplace_back(ids[a], ids[b]);
        if(res.hasWeights) res.weights.push_back(graph.weights[i]);
    }
    assert(ids.size()==maxComp.first);
    if(out) *out = ids;
    return res;
}

EdgeList undirectedToDirected(const EdgeList &graph) {
    assert(!graph.isDirected);
    auto undir = graph;
    undir.isDirected = true;
    for(int i=0; i<size(graph.edges); ++i) {
        auto [a,b] = graph.edges[i];
        undir.edges.emplace_back(b,a);
        if(graph.hasWeights) undir.weights.push_back(graph.weights[i]);
    }
    assert(size(undir.edges) == 2*size(graph.edges));
    assert(size(undir.weights) == 2*size(graph.weights));
    return undir;
}

void exportDIMACS(const EdgeList& g, const string& file, int s, int t) {
    ScopedTimer timer("export DIMACS " + file);
    ofstream f(DATA_DIR + file);
    auto m = size(g.edges);
    if(!g.isDirected) m *= 2;
    f << "p max " << g.n << ' ' << m << '\n';
    f << "n " << s+1 << " s\n";
    f << "n " << t+1 << " t\n";
    for (int i = 0; i < size(g.edges); ++i) {
        auto [a,b] = g.edges[i]; ++a; ++b;
        auto w = g.hasWeights ? g.weights[i] : 1.0;
        f << "a " << a << ' ' << b << ' ' << w << '\n';
        if(!g.isDirected)
            f << "a " << b << ' ' << a << ' ' << w << '\n';
    }
}


tuple<EdgeList,int,int> importDIMACS(const string& file) {
    ScopedTimer timer("import DIMACS " + file);
    int n,m,s=-1,t=-1;
    EdgeList g;
    g.isDirected = true;
    g.hasWeights = true;

    ifstream f(DATA_DIR + file);
    if(!f) cout << "cannot open file " << file << endl, exit(0);

    auto getType = [&] {
        string res; f>>res;
        while(res=="c") getline(f,res), f>>res;
        return res;
    };

    // read problem line
    auto cmd = getType();
    assert(cmd=="p");
    f>>cmd;
    assert(cmd=="max");
    f>>n>>m;
    g.n = n;
    g.edges.resize(m);
    g.weights.resize(m);

    // read terminal lines
    for (int i = 0; i < 2; ++i) {
        cmd = getType();
        assert(cmd=="n");
        int temp; f>>temp; --temp;
        f>>cmd;
        if(cmd=="s") s = temp;
        if(cmd=="t") t = temp;
    }
    //assert(0<=s && s<n);
    //assert(0<=t && t<n);

    // read acrs
    for (int i = 0; i < m; ++i) {
        cmd = getType();
        assert(cmd=="a");
        int a,b; f>>a>>b; --a; --b;
        assert(0<=a && a<n);
        assert(0<=b && b<n);
        g.edges[i] = {a,b};
        f >> g.weights[i];
    }

    return {g,s,t};
}

void printGraphStats(const EdgeList &graph, const string& file) {
    auto deg = degrees(graph);
    cout << "==========GRAPH============" << endl;
    if(!file.empty()) cout << "file    " << file << endl;
    cout << "n          " << graph.n << endl;
    cout << "m          " << size(graph.edges) << endl;
    cout << "weights    " << (graph.hasWeights ? "YES" : "NO") << endl;
    cout << "directed   " << (graph.isDirected ? "YES" : "NO") << endl;
    auto avg_deg = 2.0 * static_cast<double>(size(graph.edges)) / static_cast<double>(graph.n);
    auto avg_deg_o = accumulate(begin(deg),end(deg),0.0) / static_cast<double>(graph.n);
    cout << "avg deg    " <<  avg_deg << endl;
    if(avg_deg != avg_deg_o) cout << "avg deg(o) " << avg_deg_o << endl;
    cout << "min deg    " << *min_element(begin(deg),end(deg)) << endl;
    cout << "max deg    " << *max_element(begin(deg),end(deg)) << endl;
    if(!graph.isDirected) {
        auto largeCC = largestCC(graph);
        if(largeCC.n != graph.n) {
            cout << "===========LCC=============" << endl;
            cout << "n          " << largeCC.n << endl;
            cout << "m          " << size(largeCC.edges) << endl;
            cout << "avg deg    " << 2.0 * static_cast<double>(size(largeCC.edges)) / static_cast<double>(largeCC.n) << endl;
        }
    }
    cout << "===========================" << endl;
}

vector<pair<int,int>> generatePairs(const EdgeList& graph, int num, int mndeg, int mxdeg, mt19937& gen) {
    auto deg = degrees(graph);
    vector<int> candidates;
    for (int i = 0; i < size(deg); ++i) {
        if(mndeg<=deg[i] && deg[i]<=mxdeg)
            candidates.push_back(i);
    }
    if(num > size(candidates)*(size(candidates)-1)/10) {
        cout << "to few nodes in degree range!";
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
    vector res(begin(pairs), end(pairs));
    shuffle(begin(res), end(res), gen);
    return res;
}
