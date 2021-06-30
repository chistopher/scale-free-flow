
#include <BoykovKolmogorov.h>

#include "bk/maxflow.cpp"
#include "bk/graph.cpp"

using namespace std;
using namespace BoykovKolmogorov;

BKAlgorithm::BKAlgorithm(const EdgeList& graph) {
    m_graph = new Graph<double, double, double>(graph.n, static_cast<int>(size(graph.edges)));
    m_graph->add_node(graph.n);
    for (int i = 0; i < size(graph.edges); ++i) {
        auto [a, b] = graph.edges[i];
        auto w1 = graph.hasWeights ? graph.weights[i] : 1.0;
        auto w2 = graph.isDirected ? 0.0 : w1;
        m_graph->add_edge(a, b, w1, w2);
    }
}

BKAlgorithm::~BKAlgorithm() {
    delete m_graph;
}

double BKAlgorithm::maxFlow(int s, int t) {
    auto INF = 1e9;
    m_graph->add_tweights(s, INF, 0);
    m_graph->add_tweights(t, 0, INF);
    if(!first) {
        m_graph->mark_node(s);
        m_graph->mark_node(t);
    }
    auto f = m_graph->maxflow(!first); // reuse trees after first one
    m_graph->add_tweights(s, -INF, 0);
    m_graph->add_tweights(t, 0, -INF);
    m_graph->mark_node(s);
    m_graph->mark_node(t);
    return f;
}
