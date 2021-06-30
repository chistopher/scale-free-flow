
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

int main() {
    ios::sync_with_stdio(false);

    //const string file = "soc-slashdot.txt";
    //const string file = "as-skitter.txt";
    //const string file = "girg10000.txt";
    //const string file = "girglowdeg.txt";
    //const string file = "/home/chris/repos/WHFC/test_hypergraphs/Andrews.mtx.hgr.snapshot147.txt";
    const string file = "/home/chris/repos/WHFC/test_hypergraphs/Snapshot-hugebubbles-00010.mtx.hgr-2-0.txt";

    auto graph = readGraph(file);

    {
        // print stats
        cout << "===========================" << endl;
        cout << "file    " << file << endl;
        cout << "n       " << graph.n << endl;
        cout << "m       " << size(graph.edges) << endl;
        cout << "avg deg " << 2.0 * static_cast<double>(size(graph.edges)) / graph.n << endl;
        cout << "===========================" << endl;
//        graph = largestCC(graph);
        cout << "n       " << graph.n << endl;
        cout << "m       " << size(graph.edges) << endl;
        cout << "avg deg " << 2.0 * static_cast<double>(size(graph.edges)) / graph.n << endl;
        cout << "===========================" << endl;
    }

    int s = 0;
    int t = 2*6752+1;
    t = 2*1904564+1;
    {
        ScopedTimer timer("bi");
        cout << Dinics::Dinics4Skip(graph).maxFlow(s,t) << endl;
    }
    {
        ScopedTimer timer("uni");
        cout << Dinics::Dinics0Vanilla(graph).maxFlow(s,t) << endl;
    }

    return 0;
}
